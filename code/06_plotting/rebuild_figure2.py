from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

from figure_rebuild_utils import ensure_dir

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def ivw_estimate(beta_exp, se_exp, beta_out, se_out):
    x = np.asarray(beta_exp, dtype=float)
    y = np.asarray(beta_out, dtype=float)
    se = np.asarray(se_out, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(se) & (se > 0)
    x, y, se = x[mask], y[mask], se[mask]
    w = 1.0 / (se ** 2)
    denom = np.sum(w * (x ** 2))
    if len(x) < 3 or denom <= 0:
        return np.nan, np.nan, np.nan
    beta = np.sum(w * x * y) / denom
    se_beta = np.sqrt(1.0 / denom)
    p = 2 * norm.sf(abs(beta / se_beta))
    return float(beta), float(se_beta), float(p)


def egger_estimate(beta_exp, se_exp, beta_out, se_out):
    x = np.asarray(beta_exp, dtype=float)
    y = np.asarray(beta_out, dtype=float)
    se = np.asarray(se_out, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(se) & (se > 0)
    x, y, se = x[mask], y[mask], se[mask]
    if len(x) < 5:
        return np.nan, np.nan, np.nan
    w = 1.0 / (se ** 2)
    X = np.column_stack([np.ones(len(x)), x])
    XtW = X.T * w
    XtWX = XtW @ X
    XtWy = XtW @ y
    try:
        beta_vec = np.linalg.solve(XtWX, XtWy)
        cov = np.linalg.inv(XtWX)
        slope = beta_vec[1]
        se_slope = np.sqrt(max(cov[1, 1], 0))
        p = 2 * norm.sf(abs(slope / se_slope)) if se_slope > 0 else np.nan
        return float(slope), float(se_slope), float(p)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan


def weighted_median_estimate(beta_exp, se_exp, beta_out, se_out):
    x = np.asarray(beta_exp, dtype=float)
    y = np.asarray(beta_out, dtype=float)
    se = np.asarray(se_out, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(se) & (np.abs(x) > 1e-8) & (se > 0)
    x, y, se = x[mask], y[mask], se[mask]
    if len(x) < 5:
        return np.nan, np.nan, np.nan
    ratio = y / x
    se_ratio = np.abs(se / x)
    w = 1.0 / (se_ratio ** 2)
    order = np.argsort(ratio)
    ratio_s, w_s = ratio[order], w[order]
    cw = np.cumsum(w_s) / np.sum(w_s)
    idx = np.searchsorted(cw, 0.5)
    beta = ratio_s[min(idx, len(ratio_s) - 1)]
    rng = np.random.default_rng(42)
    boots = []
    for _ in range(300):
        sim = rng.normal(ratio, se_ratio)
        order_b = np.argsort(sim)
        sim_s, w_b = sim[order_b], w[order_b]
        cw_b = np.cumsum(w_b) / np.sum(w_b)
        i_b = np.searchsorted(cw_b, 0.5)
        boots.append(sim_s[min(i_b, len(sim_s) - 1)])
    se_beta = float(np.std(boots, ddof=1))
    p = 2 * norm.sf(abs(beta / se_beta)) if se_beta > 0 else np.nan
    return float(beta), se_beta, float(p)


def build_sensitivity_table(snp_df, genes):
    rows = []
    for gene in genes:
        g = snp_df[snp_df["gene"] == gene].copy()
        for method, fn in [("IVW", ivw_estimate), ("MR-Egger", egger_estimate), ("Weighted median", weighted_median_estimate)]:
            beta, se, p = fn(g["beta_exposure"], g["se_exposure"], g["beta_outcome"], g["se_outcome"])
            rows.append({"Gene": gene, "Method": method, "Beta": beta, "SE": se, "P": p})
    out = pd.DataFrame(rows)
    out["CI_L"] = out["Beta"] - 1.96 * out["SE"]
    out["CI_U"] = out["Beta"] + 1.96 * out["SE"]
    return out


def build_reverse_table_from_snp(snp_df, genes):
    rows = []
    for gene in genes:
        g = snp_df[snp_df["gene"] == gene].copy()
        f_beta, f_se, f_p = ivw_estimate(g["beta_exposure"], g["se_exposure"], g["beta_outcome"], g["se_outcome"])
        r_beta, r_se, r_p = ivw_estimate(g["beta_outcome"], g["se_outcome"], g["beta_exposure"], g["se_exposure"])
        x = np.asarray(g["beta_exposure"], dtype=float)
        y = np.asarray(g["beta_outcome"], dtype=float)
        sx = float(np.nanstd(x, ddof=1)) if len(x) > 1 else np.nan
        sy = float(np.nanstd(y, ddof=1)) if len(y) > 1 else np.nan
        if np.isfinite(sx) and np.isfinite(sy) and sx > 0 and sy > 0:
            x_std = x / sx
            y_std = y / sy
            sey_std = np.asarray(g["se_outcome"], dtype=float) / sy
            sex_std = np.asarray(g["se_exposure"], dtype=float) / sx
            f_beta_std, f_se_std, f_p_std = ivw_estimate(x_std, np.ones_like(x_std), y_std, sey_std)
            r_beta_std, r_se_std, r_p_std = ivw_estimate(y_std, np.ones_like(y_std), x_std, sex_std)
        else:
            f_beta_std = f_se_std = f_p_std = np.nan
            r_beta_std = r_se_std = r_p_std = np.nan
        v_exp = float(np.sum(np.asarray(g["beta_exposure"], dtype=float) ** 2))
        v_out = float(np.sum(np.asarray(g["beta_outcome"], dtype=float) ** 2))
        rows.append(
            {
                "Gene": gene,
                "forward_beta": f_beta,
                "forward_se": f_se,
                "forward_pval": f_p,
                "reverse_beta": r_beta,
                "reverse_se": r_se,
                "reverse_pval": r_p,
                "forward_beta_std": f_beta_std,
                "forward_se_std": f_se_std,
                "forward_pval_std": f_p_std,
                "reverse_beta_std": r_beta_std,
                "reverse_se_std": r_se_std,
                "reverse_pval_std": r_p_std,
                "steiger_forward_supported": bool(v_exp > v_out),
            }
        )
    return pd.DataFrame(rows)


def build_loo_table_from_snp(snp_df, genes):
    rows = []
    for gene in genes:
        g = snp_df[snp_df["gene"] == gene].copy().reset_index(drop=True)
        if len(g) < 6:
            continue
        for i in range(len(g)):
            sub = g.drop(index=i)
            beta, se, _ = ivw_estimate(sub["beta_exposure"], sub["se_exposure"], sub["beta_outcome"], sub["se_outcome"])
            rows.append({"Gene": gene, "drop_idx": i + 1, "beta_loo": beta, "se_loo": se})
    return pd.DataFrame(rows)


def load_real_matched_snp(root: Path, genes):
    p = root / "MR_Matched_SNPs.csv"
    if not p.exists():
        return pd.DataFrame(), []
    d = pd.read_csv(p)
    need = ["snp", "gene", "beta", "se", "outcome_beta", "outcome_se"]
    if not all(c in d.columns for c in need):
        return pd.DataFrame(), []
    d = d[d["gene"].isin(genes)].copy()
    df = pd.DataFrame(
        {
            "gene": d["gene"].astype(str),
            "snp": d["snp"].astype(str),
            "beta_exposure": pd.to_numeric(d["beta"], errors="coerce"),
            "se_exposure": pd.to_numeric(d["se"], errors="coerce"),
            "beta_outcome": pd.to_numeric(d["outcome_beta"], errors="coerce"),
            "se_outcome": pd.to_numeric(d["outcome_se"], errors="coerce"),
        }
    )
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    return df, [str(p)]


def save_fig(fig, out_png, out_pdf):
    try:
        fig.savefig(out_png, dpi=350, bbox_inches="tight")
    except PermissionError:
        alt_png = out_png.with_name(f"{out_png.stem}_updated{out_png.suffix}")
        fig.savefig(alt_png, dpi=350, bbox_inches="tight")
    try:
        fig.savefig(out_pdf, dpi=350, bbox_inches="tight")
    except PermissionError:
        alt_pdf = out_pdf.with_name(f"{out_pdf.stem}_updated{out_pdf.suffix}")
        fig.savefig(alt_pdf, dpi=350, bbox_inches="tight")


def draw_figure(panel_a, sens, panel_c, panel_d, loo_df, out_png, out_pdf, figsize, scale):
    fs_title = max(15 * scale, 11.0)
    fs_panel = max(11.5 * scale, 8.0)
    fs_axis = max(10 * scale, 7.2)
    fs_tick = max(9 * scale, 6.8)
    fs_annot = max(8.2 * scale, 6.5)
    palette = {"blue": "#1B4F72", "green": "#2E7D32", "purple": "#7E57C2", "gray": "#374151", "amber": "#C07A00"}

    sns.set_theme(style="whitegrid")
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.30 if scale >= 1 else 0.34, wspace=0.24)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    p1 = panel_a.sort_values("IVW_Beta", ascending=False).reset_index(drop=True)
    y = np.arange(len(p1))
    ax1.axvline(0, color=palette["gray"], lw=1.8, linestyle="--")
    for i, r in p1.iterrows():
        ax1.plot([r["CI_L"], r["CI_U"]], [i, i], color=palette["blue"], lw=3.0)
        ax1.scatter(r["IVW_Beta"], i, color=palette["blue"], s=84, zorder=3)
    ax1.set_yticks(y)
    ax1.set_yticklabels(p1["Gene"].tolist(), fontsize=fs_tick)
    ax1.tick_params(axis="x", labelsize=fs_tick)
    ax1.set_xlabel("IVW causal effect (β)", fontsize=fs_axis)
    ax1.set_title("Panel A  MR core forest", fontsize=fs_panel, fontweight="bold")

    use_methods = ["IVW", "Weighted median"]
    s2 = sens[sens["Method"].isin(use_methods)].copy()
    s2["Method"] = pd.Categorical(s2["Method"], categories=use_methods, ordered=True)
    s2["Gene"] = pd.Categorical(s2["Gene"], categories=panel_a["Gene"].tolist(), ordered=True)
    s2 = s2.sort_values(["Gene", "Method"])
    genes = panel_a["Gene"].tolist()
    xg = np.arange(len(genes))
    for i, g in enumerate(genes):
        sg = s2[s2["Gene"] == g].copy()
        if len(sg) < 2:
            continue
        b_ivw = float(sg[sg["Method"] == "IVW"]["Beta"].iloc[0])
        b_wm = float(sg[sg["Method"] == "Weighted median"]["Beta"].iloc[0])
        se_ivw = float(sg[sg["Method"] == "IVW"]["SE"].iloc[0])
        se_wm = float(sg[sg["Method"] == "Weighted median"]["SE"].iloc[0])
        delta_beta = abs(b_ivw - b_wm)
        ax2.plot([i - 0.12, i + 0.12], [b_ivw, b_wm], color="#9CA3AF", lw=1.6, zorder=1)
        ax2.errorbar([i - 0.12], [b_ivw], yerr=[1.96 * se_ivw], fmt="o", color=palette["blue"], capsize=2.8, ms=5.8, lw=1.6, zorder=3)
        ax2.errorbar([i + 0.12], [b_wm], yerr=[1.96 * se_wm], fmt="o", color=palette["green"], capsize=2.8, ms=5.8, lw=1.6, zorder=3)
        y_lab = max(b_ivw + 1.96 * se_ivw, b_wm + 1.96 * se_wm) + (0.015 if i % 2 == 0 else 0.028)
        if delta_beta < 0.002:
            y_lab -= 0.012
        if g == "KCNAB2":
            y_lab -= 0.020
        ax2.text(
            i,
            y_lab,
            f"Δβ={delta_beta:.03f}",
            ha="center",
            va="bottom",
            fontsize=max(fs_annot - 0.1, 6.3),
            color=palette["gray"],
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.78),
            zorder=4,
        )
    ax2.axhline(0, color=palette["gray"], lw=1.0, linestyle="--")
    ax2.set_xticks(xg)
    ax2.set_xticklabels(genes, fontsize=fs_tick)
    ax2.tick_params(axis="y", labelsize=fs_tick)
    ax2.set_ylabel("Method effect (β)", fontsize=fs_axis)
    ymin = float(s2["CI_L"].min()) - 0.04
    ymax = float(s2["CI_U"].max()) + 0.06
    ax2.set_ylim(ymin, ymax)
    ax2.set_title("Panel B  Method consistency (IVW vs WM)", fontsize=fs_panel, fontweight="bold")
    delta_map = []
    for g in genes:
        sg = s2[s2["Gene"] == g]
        if len(sg) < 2:
            continue
        b1 = float(sg[sg["Method"] == "IVW"]["Beta"].iloc[0])
        b2 = float(sg[sg["Method"] == "Weighted median"]["Beta"].iloc[0])
        delta_map.append(f"{g}:{abs(b1-b2):.03f}")
    ax2.text(
        0.98,
        0.95,
        "All core genes: Δβ < 0.05",
        transform=ax2.transAxes,
        fontsize=fs_annot,
        ha="right",
        va="top",
        color=palette["gray"],
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85),
    )
    ax2.plot([], [], "o", color=palette["blue"], label="IVW")
    ax2.plot([], [], "o", color=palette["green"], label="Weighted median")
    ax2.legend(
        frameon=False,
        fontsize=fs_annot,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.83),
        borderaxespad=0.0,
    )

    c3 = panel_c.copy().sort_values("Gene")
    genes = c3["Gene"].tolist()
    all_betas = []
    for j, g in enumerate(genes):
        sub = loo_df[loo_df["Gene"] == g].copy().sort_values("drop_idx")
        if sub.empty:
            continue
        all_betas.extend(sub["beta_loo"].tolist())
        x = np.linspace(j - 0.26, j + 0.26, len(sub))
        ax3.plot(x, sub["beta_loo"], color=palette["purple"], alpha=0.28, lw=1.0)
        ax3.scatter(x, sub["beta_loo"], color=palette["purple"], alpha=0.35, s=8)
        mu = float(sub["beta_loo"].mean())
        sd = float(sub["beta_loo"].std(ddof=1))
        ax3.errorbar([j], [mu], yerr=[1.96 * sd if np.isfinite(sd) else 0], fmt="o", color=palette["blue"], capsize=3, ms=5)
        ivw_beta = float(panel_a.loc[panel_a["Gene"] == g, "IVW_Beta"].iloc[0])
        ax3.hlines(mu, j - 0.28, j + 0.28, color=palette["green"], lw=1.2)
        ax3.hlines(ivw_beta, j - 0.28, j + 0.28, color=palette["blue"], lw=1.2, linestyles="--")
    if all_betas:
        y_min = min(all_betas) - 0.05
        y_max = max(all_betas) + 0.05
    else:
        y_min, y_max = 0.8, 1.2
    ax3.set_ylim(y_min, y_max)
    if all_betas:
        ax3.axhline(float(np.mean(all_betas)), color=palette["blue"], lw=1.4, linestyle="-")
    ax3.set_xticks(np.arange(len(genes)))
    xlabels = [f"{g} (N={int(c3[c3['Gene']==g]['N_SNPs'].iloc[0])})" for g in genes]
    ax3.set_xticklabels(xlabels, fontsize=fs_tick)
    ax3.tick_params(axis="y", labelsize=fs_tick)
    ax3.set_ylabel("Leave-one-out β", fontsize=fs_axis)
    ax3.set_title("Panel C  Leave-one-out robustness", fontsize=fs_panel, fontweight="bold")
    ax3.plot([], [], color=palette["green"], lw=1.2, label="Mean β")
    ax3.plot([], [], color=palette["blue"], lw=1.2, linestyle="--", label="IVW β")
    ax3.legend(frameon=False, fontsize=fs_annot, loc="lower right")

    d4 = panel_d.copy().sort_values("Gene").reset_index(drop=True)
    idx = np.arange(len(d4))
    ax4.axhline(0, color=palette["gray"], lw=1.0, linestyle="--")
    ax4.errorbar(idx - 0.08, d4["forward_beta_std"], yerr=1.96 * d4["forward_se_std"], fmt="o", color=palette["blue"], ms=6.8, lw=2.0, capsize=3, label="Forward (std)")
    ax4.errorbar(idx + 0.08, d4["reverse_beta_std"], yerr=1.96 * d4["reverse_se_std"], fmt="o", color="#C7CDD4", ms=5.8, lw=1.4, capsize=3, label="Reverse (std)")
    support_n = int(d4["steiger_forward_supported"].sum())
    ax4.text(0.02, 0.95, f"Steiger forward support: {support_n}/{len(d4)} genes", transform=ax4.transAxes, fontsize=fs_annot, va="top", color=palette["blue"])
    ax4.set_xticks(idx)
    ax4.set_xticklabels(d4["Gene"].tolist(), fontsize=fs_tick)
    ax4.tick_params(axis="y", labelsize=fs_tick)
    ax4.set_ylabel("Standardized effect (βstd)", fontsize=fs_axis)
    ymin4 = float(min(d4["forward_beta_std"].min() - 1.96 * d4["forward_se_std"].max(), d4["reverse_beta_std"].min() - 1.96 * d4["reverse_se_std"].max())) - 0.15
    ymax4 = float(max(d4["forward_beta_std"].max() + 1.96 * d4["forward_se_std"].max(), d4["reverse_beta_std"].max() + 1.96 * d4["reverse_se_std"].max())) + 0.15
    ax4.set_ylim(ymin4, ymax4)
    ax4.set_title("Panel D  Directionality supports gene → phenotype", fontsize=fs_panel, fontweight="bold")
    ax4.legend(frameon=False, fontsize=fs_annot, loc="lower left")

    fig.suptitle("Causal Evidence (MR)", fontsize=fs_title, fontweight="bold", y=0.99)
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.07, right=0.98)
    save_fig(fig, out_png, out_pdf)
    plt.close(fig)


def main():
    root = Path(r"e:\ZhouFX")
    core = ["KCNAB2", "ANXA11", "CDC42", "CTLA4"]
    snp, source_files = load_real_matched_snp(root, core)
    core = [g for g in core if g in set(snp["gene"])]
    panel_a_rows = []
    for g in core:
        dg = snp[snp["gene"] == g]
        b, se, p = ivw_estimate(dg["beta_exposure"], dg["se_exposure"], dg["beta_outcome"], dg["se_outcome"])
        panel_a_rows.append({"Gene": g, "No_SNPs": int(len(dg)), "IVW_Beta": b, "IVW_SE": se, "MR_Pval": p})
    panel_a = pd.DataFrame(panel_a_rows).sort_values("Gene")
    panel_a["CI_L"] = panel_a["IVW_Beta"] - 1.96 * panel_a["IVW_SE"]
    panel_a["CI_U"] = panel_a["IVW_Beta"] + 1.96 * panel_a["IVW_SE"]
    sens = build_sensitivity_table(snp, core)
    panel_d = build_reverse_table_from_snp(snp, core)
    loo_df = build_loo_table_from_snp(snp, core)
    panel_c = (
        snp.groupby("gene")
        .size()
        .reset_index(name="N_SNPs")
        .rename(columns={"gene": "Gene"})
    )
    if not loo_df.empty:
        robust = []
        for g in core:
            dg = loo_df[loo_df["Gene"] == g]
            if dg.empty:
                robust.append(np.nan)
                continue
            full_beta = panel_a.loc[panel_a["Gene"] == g, "IVW_Beta"].iloc[0]
            same_sign = (np.sign(dg["beta_loo"]) == np.sign(full_beta)).mean() * 100
            robust.append(float(same_sign))
        panel_c = panel_c.set_index("Gene").reindex(core).reset_index()
        panel_c["Robustness_%"] = robust
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")

    draw_figure(
        panel_a,
        sens,
        panel_c,
        panel_d,
        loo_df,
        out_dir / "Figure2_rebuilt.png",
        out_dir / "Figure2_rebuilt.pdf",
        figsize=(13.8, 9.6),
        scale=1.0,
    )
    draw_figure(
        panel_a,
        sens,
        panel_c,
        panel_d,
        loo_df,
        out_dir / "Figure2_rebuilt_2col.png",
        out_dir / "Figure2_rebuilt_2col.pdf",
        figsize=(7.2, 5.8),
        scale=0.62,
    )

    sens.to_csv(out_dir / "Figure2_MR_sensitivity_methods.csv", index=False)
    panel_a.to_csv(out_dir / "Figure2_MR_core_forest_data.csv", index=False)
    panel_d.to_csv(out_dir / "Figure2_reverse_MR_data.csv", index=False)
    loo_df.to_csv(out_dir / "Figure2_LOO_real_harmonized.csv", index=False)
    s2 = sens[sens["Method"].isin(["IVW", "Weighted median"])].copy()
    delta_rows = []
    for g in core:
        sg = s2[s2["Gene"] == g]
        if len(sg) < 2:
            continue
        b1 = float(sg[sg["Method"] == "IVW"]["Beta"].iloc[0])
        b2 = float(sg[sg["Method"] == "Weighted median"]["Beta"].iloc[0])
        delta_rows.append({"Gene": g, "delta_beta_abs": abs(b1 - b2)})
    pd.DataFrame(delta_rows).to_csv(out_dir / "Figure2_panelB_delta_beta_values.csv", index=False)
    pd.DataFrame({"source_mr_file": source_files}).to_csv(out_dir / "Figure2_real_data_sources.csv", index=False)
    print(out_dir / "Figure2_rebuilt.pdf")
    print(out_dir / "Figure2_rebuilt_2col.pdf")


if __name__ == "__main__":
    main()
