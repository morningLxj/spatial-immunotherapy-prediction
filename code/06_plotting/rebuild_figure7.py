from pathlib import Path
import gzip
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from figure_rebuild_utils import ensure_dir

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def safe_savefig(fig, path, **kwargs):
    try:
        fig.savefig(path, **kwargs)
        return path
    except PermissionError:
        alt = path.with_name(f"{path.stem}_updated{path.suffix}")
        fig.savefig(alt, **kwargs)
        return alt


def load_top_features(root: Path, n=120):
    p = root / "multi_omics_feature_importance.csv"
    feat = pd.read_csv(p)
    features = feat["feature"].head(n).tolist()
    ens_base = [f.split("_")[0].split(".")[0] for f in features]
    return features, ens_base


def load_tcga(root: Path):
    df = pd.read_csv(root / "integrated_data.csv")
    req = ["sample_id", "OS.time", "OS", "response"]
    for c in req:
        if c not in df.columns:
            raise ValueError(f"Missing {c}")
    df = df.dropna(subset=["OS.time", "OS", "response"]).copy()
    df = df[df["OS.time"] > 0].copy()
    df["OS"] = df["OS"].astype(int)
    df["response"] = df["response"].astype(int)
    return df


def load_gse31210(root: Path, top_features, ens_base):
    clinical = pd.read_csv(root / "gse31210_clinical.csv")
    if "GSM_ID" not in clinical.columns:
        raise ValueError("GSM_ID missing in gse31210_clinical.csv")
    clinical = clinical.set_index("GSM_ID")
    exc = "exclude for prognosis analysis due to incomplete resection or adjuvant therapy"
    if exc in clinical.columns:
        clinical = clinical[clinical[exc] != "exclude"]
    if "gene alteration status" in clinical.columns:
        clinical = clinical[~clinical["gene alteration status"].isin(["KRAS mutation +", "ALK-fusion +"])]
    clinical = clinical.dropna(subset=["OS.time", "OS"]).copy()
    clinical["OS.time"] = pd.to_numeric(clinical["OS.time"], errors="coerce")
    clinical["OS"] = pd.to_numeric(clinical["OS"], errors="coerce")
    clinical = clinical.dropna(subset=["OS.time", "OS"])
    clinical = clinical[clinical["OS.time"] > 0]
    clinical["OS"] = clinical["OS"].astype(int)

    map_df = pd.read_csv(root / "gpl570_probe_mapping.csv")
    map_df = map_df.dropna(subset=["Probe_ID", "Ensembl_ID"]).copy()
    map_df["Ensembl_ID"] = map_df["Ensembl_ID"].astype(str).str.split(".").str[0]
    target = set(ens_base)
    map_df = map_df[map_df["Ensembl_ID"].isin(target)].copy()
    probe_to_ens = dict(zip(map_df["Probe_ID"].astype(str), map_df["Ensembl_ID"]))
    target_probes = set(probe_to_ens.keys())

    expr_file = root / "GSE31210_series_matrix.txt.gz"
    rows = []
    cols = None
    with gzip.open(expr_file, "rt", encoding="utf-8") as f:
        in_table = False
        for line in f:
            if line.strip() == "!series_matrix_table_begin":
                in_table = True
                header_line = next(f)
                cols = [x.strip('"') for x in header_line.strip().split("\t")]
                continue
            if in_table:
                if line.strip() == "!series_matrix_table_end":
                    break
                parts = line.strip().split("\t")
                probe = parts[0].strip('"')
                if probe in target_probes:
                    vals = [float(x.strip('"')) if x.strip('"') else np.nan for x in parts[1:]]
                    rows.append([probe] + vals)
    if cols is None or len(rows) == 0:
        raise ValueError("Failed to parse GSE31210 matrix")
    expr = pd.DataFrame(rows, columns=["ID_REF"] + cols[1:]).set_index("ID_REF")
    if np.nanmax(expr.values) > 20:
        expr = np.log2(expr + 1)
    expr["Ensembl_ID"] = expr.index.map(probe_to_ens)
    gene_expr = expr.groupby("Ensembl_ID").mean().T
    ens_to_feat = {}
    for f in top_features:
        ens_to_feat[f.split("_")[0].split(".")[0]] = f
    gene_expr = gene_expr.rename(columns={c: ens_to_feat.get(c, c) for c in gene_expr.columns})
    common_samples = gene_expr.index.intersection(clinical.index)
    out = clinical.loc[common_samples].copy().join(gene_expr.loc[common_samples], how="left")
    out = out.dropna(subset=["OS.time", "OS"]).copy()
    return out


def train_predict_scores(tcga_df, gse_df, features):
    common = [f for f in features if (f in tcga_df.columns) and (f in gse_df.columns)]
    if len(common) < 20:
        raise ValueError("Too few shared features for cross-platform validation")
    X_tcga = tcga_df[common]
    y_tcga = tcga_df["response"].astype(int).values
    imp = SimpleImputer(strategy="median")
    X_tcga_imp = imp.fit_transform(X_tcga)
    sc = StandardScaler()
    X_tcga_sc = sc.fit_transform(X_tcga_imp)
    clf = XGBClassifier(
        random_state=42,
        n_estimators=200,
        learning_rate=0.05,
        max_depth=4,
        subsample=0.9,
        colsample_bytree=0.9,
        eval_metric="logloss",
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    tcga_scores = cross_val_predict(clf, X_tcga_sc, y_tcga, cv=cv, method="predict_proba")[:, 1]
    clf.fit(X_tcga_sc, y_tcga)
    X_gse = gse_df[common]
    X_gse_imp = imp.transform(X_gse)
    X_gse_sc = sc.transform(X_gse_imp)
    gse_scores = clf.predict_proba(X_gse_sc)[:, 1]
    return tcga_scores, gse_scores


def km_and_hr(df, score_col, time_col, event_col):
    d = df.copy()
    cutoff = float(np.median(d[score_col].values))
    d["Risk_Group"] = np.where(d[score_col] > cutoff, "High Risk", "Low Risk")
    t = d[time_col].values.astype(float) / 365.0
    e = d[event_col].values.astype(int)
    high = d["Risk_Group"].values == "High Risk"
    low = ~high
    lr = logrank_test(t[high], t[low], event_observed_A=e[high], event_observed_B=e[low])
    cph = CoxPHFitter()
    cph_df = pd.DataFrame({"T": t, "E": e, "G": high.astype(int)})
    cph.fit(cph_df, duration_col="T", event_col="E")
    hr = float(cph.summary.loc["G", "exp(coef)"])
    l95 = float(cph.summary.loc["G", "exp(coef) lower 95%"])
    u95 = float(cph.summary.loc["G", "exp(coef) upper 95%"])
    se = float(cph.summary.loc["G", "se(coef)"])
    return d, cutoff, float(lr.p_value), hr, l95, u95, se


def pooled_meta(stats_df):
    yi = np.log(stats_df["hr"].values)
    vi = stats_df["se"].values ** 2
    wi = 1.0 / vi
    mu_fix = float(np.sum(wi * yi) / np.sum(wi))
    q = float(np.sum(wi * (yi - mu_fix) ** 2))
    k = len(yi)
    c = float(np.sum(wi) - np.sum(wi**2) / np.sum(wi))
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 else 0.0
    wi_re = 1.0 / (vi + tau2)
    mu_re = float(np.sum(wi_re * yi) / np.sum(wi_re))
    se_re = float(np.sqrt(1.0 / np.sum(wi_re)))
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > (k - 1) and q > 0 else 0.0
    return np.exp(mu_re), np.exp(mu_re - 1.96 * se_re), np.exp(mu_re + 1.96 * se_re), i2


def _binary_target_at_time(df, years, time_col, event_col):
    t_days = years * 365.0
    valid = ~((df[event_col].values == 0) & (df[time_col].values <= t_days))
    y = ((df.loc[valid, event_col].values == 1) & (df.loc[valid, time_col].values <= t_days)).astype(int)
    return valid, y


def auc_with_ci(df, scores, years, time_col, event_col, n_boot=300):
    valid, y = _binary_target_at_time(df, years, time_col, event_col)
    s = np.asarray(scores)[valid]
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan
    fpr, tpr, _ = roc_curve(y, s)
    auc0 = float(auc(fpr, tpr))
    rng = np.random.default_rng(42 + years)
    boots = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            continue
        sb = s[idx]
        fprb, tprb, _ = roc_curve(yb, sb)
        boots.append(float(auc(fprb, tprb)))
    if len(boots) == 0:
        return auc0, np.nan, np.nan
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return auc0, float(lo), float(hi)


def decile_event_rates(df, scores, time_col, event_col):
    d = df.copy()
    d["score"] = np.asarray(scores)
    d["decile"] = pd.qcut(d["score"], 10, labels=False, duplicates="drop")
    g = d.groupby("decile", as_index=False).agg(event_rate=(event_col, "mean"), score_mean=("score", "mean"))
    return g


def draw_figure7(tcga_df, tcga_scores, gse_df, gse_scores, out_png, out_pdf, figsize, scale):
    fs_panel = max(11.0 * scale, 8.2)
    fs_axis = max(9.8 * scale, 7.2)
    fs_tick = max(8.4 * scale, 6.4)
    fs_note = max(7.8 * scale, 6.0)
    fs_title = max(14.8 * scale, 11.0)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.30)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    gse_plot = gse_df.copy()
    gse_plot["score"] = gse_scores
    d, _, p, hr, l95, u95, _ = km_and_hr(gse_plot, "score", "OS.time", "OS")
    t = d["OS.time"].values / 365.0
    e = d["OS"].values.astype(int)
    high = d["Risk_Group"] == "High Risk"
    low = d["Risk_Group"] == "Low Risk"
    kmh = KaplanMeierFitter()
    kml = KaplanMeierFitter()
    kmh.fit(t[high], event_observed=e[high], label="High risk")
    kml.fit(t[low], event_observed=e[low], label="Low risk")
    kmh.plot_survival_function(ax=axA, ci_show=False, color="#D55E00", linewidth=2.0)
    kml.plot_survival_function(ax=axA, ci_show=False, color="#3C5488", linewidth=2.0)
    axA.set_title("Representative external cohort survival separation", fontsize=fs_panel, fontweight="bold")
    axA.set_xlabel("Time (years)", fontsize=fs_axis)
    axA.set_ylabel("Survival probability", fontsize=fs_axis)
    axA.tick_params(labelsize=fs_tick)
    axA.legend(frameon=False, fontsize=fs_note, loc="upper right")
    axA.text(
        0.03,
        0.95,
        f"GSE31210: HR={hr:.2f} (95% CI {l95:.2f}-{u95:.2f})\nlog-rank p={p:.3g}",
        transform=axA.transAxes,
        ha="left",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )

    tcga_plot = tcga_df.copy()
    tcga_plot["score"] = tcga_scores
    cohorts = []
    _, _, p1, hr1, l1, u1, se1 = km_and_hr(tcga_plot, "score", "OS.time", "OS")
    cohorts.append({"cohort": "TCGA (internal)", "hr": hr1, "lcl": l1, "ucl": u1, "p": p1, "se": se1})
    _, _, p2, hr2, l2, u2, se2 = km_and_hr(gse_plot, "score", "OS.time", "OS")
    cohorts.append({"cohort": "GSE31210 (external)", "hr": hr2, "lcl": l2, "ucl": u2, "p": p2, "se": se2})
    if "pstage iorii" in gse_plot.columns:
        stage_col = gse_plot["pstage iorii"].astype(str)
        m1 = stage_col.str.contains("I", regex=False) & (~stage_col.str.contains("II", regex=False))
        m2 = stage_col.str.contains("II", regex=False)
        if m1.sum() >= 30:
            d1 = gse_plot.loc[m1].copy()
            _, _, p3, hr3, l3, u3, se3 = km_and_hr(d1, "score", "OS.time", "OS")
            cohorts.append({"cohort": "GSE31210 Stage I", "hr": hr3, "lcl": l3, "ucl": u3, "p": p3, "se": se3})
        if m2.sum() >= 20:
            d2 = gse_plot.loc[m2].copy()
            _, _, p4, hr4, l4, u4, se4 = km_and_hr(d2, "score", "OS.time", "OS")
            cohorts.append({"cohort": "GSE31210 Stage II", "hr": hr4, "lcl": l4, "ucl": u4, "p": p4, "se": se4})
    forest = pd.DataFrame(cohorts)
    pooled_hr, pooled_l, pooled_u, i2 = pooled_meta(forest)
    forest_plot = pd.concat(
        [forest, pd.DataFrame([{"cohort": "Pooled (random-effects)", "hr": pooled_hr, "lcl": pooled_l, "ucl": pooled_u, "p": np.nan, "se": np.nan}])],
        ignore_index=True,
    )
    y = np.arange(len(forest_plot))[::-1]
    axB.hlines(y, forest_plot["lcl"], forest_plot["ucl"], color="#3C5488", lw=1.8)
    axB.scatter(forest_plot["hr"], y, color="#D55E00", s=30, zorder=3)
    axB.axvline(1.0, color="black", lw=1, ls="--")
    axB.set_xscale("log")
    axB.set_yticks(y)
    wrapped_labels = []
    for lbl in forest_plot["cohort"].tolist():
        if lbl == "TCGA (internal)":
            wrapped_labels.append("TCGA\n(internal)")
        elif lbl == "GSE31210 (external)":
            wrapped_labels.append("GSE31210\n(external)")
        elif lbl == "GSE31210 Stage I":
            wrapped_labels.append("GSE31210\nStage I")
        elif lbl == "GSE31210 Stage II":
            wrapped_labels.append("GSE31210\nStage II")
        elif lbl == "Pooled (random-effects)":
            wrapped_labels.append("Pooled\n(random-effects)")
        else:
            wrapped_labels.append(lbl)
    axB.set_yticklabels(wrapped_labels, fontsize=fs_tick, linespacing=1.05)
    axB.set_xlabel("Hazard ratio (log scale)", fontsize=fs_axis)
    axB.set_title("Cross-cohort hazard consistency with pooled effect", fontsize=fs_panel, fontweight="bold")
    axB.tick_params(labelsize=fs_tick)
    axB.text(
        0.97,
        0.95,
        f"Pooled HR={pooled_hr:.2f} ({pooled_l:.2f}-{pooled_u:.2f})\nI²={i2:.1f}%",
        transform=axB.transAxes,
        ha="right",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )

    tcga_dec = decile_event_rates(tcga_plot, tcga_scores, "OS.time", "OS")
    gse_dec = decile_event_rates(gse_plot, gse_scores, "OS.time", "OS")
    k = int(min(len(tcga_dec), len(gse_dec)))
    x = tcga_dec["event_rate"].values[:k]
    yv = gse_dec["event_rate"].values[:k]
    axC.scatter(x, yv, s=30, color="#009E73", alpha=0.9)
    lim_lo = float(min(np.min(x), np.min(yv)) * 0.9)
    lim_hi = float(max(np.max(x), np.max(yv)) * 1.1 + 1e-6)
    axC.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", lw=1)
    rho, pval = stats.spearmanr(x, yv)
    axC.set_xlim(lim_lo, lim_hi)
    axC.set_ylim(lim_lo, lim_hi)
    axC.set_xlabel("TCGA decile event rate", fontsize=fs_axis)
    axC.set_ylabel("GSE31210 decile event rate", fontsize=fs_axis)
    axC.set_title("Cross-platform risk-gradient consistency", fontsize=fs_panel, fontweight="bold")
    axC.tick_params(labelsize=fs_tick)
    axC.text(
        0.97,
        0.95,
        f"Spearman ρ={rho:.2f}, p={pval:.3g}",
        transform=axC.transAxes,
        ha="right",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )

    rows = []
    for cohort_name, dfx, sx in [("TCGA", tcga_plot, tcga_scores), ("GSE31210", gse_plot, gse_scores)]:
        for yr in [1, 3, 5]:
            a0, lo, hi = auc_with_ci(dfx, sx, yr, "OS.time", "OS", n_boot=250)
            rows.append({"cohort": cohort_name, "year": yr, "auc": a0, "lo": lo, "hi": hi})
    met = pd.DataFrame(rows)
    x_positions = np.arange(3)
    width = 0.34
    tc = met[met["cohort"] == "TCGA"].sort_values("year")
    ge = met[met["cohort"] == "GSE31210"].sort_values("year")
    axD.bar(x_positions - width / 2, tc["auc"].values, width=width, color="#3C5488", alpha=0.85, label="TCGA")
    axD.bar(x_positions + width / 2, ge["auc"].values, width=width, color="#D55E00", alpha=0.85, label="GSE31210")
    axD.errorbar(x_positions - width / 2, tc["auc"].values, yerr=[tc["auc"].values - tc["lo"].values, tc["hi"].values - tc["auc"].values], fmt="none", ecolor="black", capsize=2, lw=1)
    axD.errorbar(x_positions + width / 2, ge["auc"].values, yerr=[ge["auc"].values - ge["lo"].values, ge["hi"].values - ge["auc"].values], fmt="none", ecolor="black", capsize=2, lw=1)
    axD.set_xticks(x_positions)
    axD.set_xticklabels(["1-year", "3-year", "5-year"], fontsize=fs_tick)
    axD.set_ylim(0.45, 0.82)
    axD.set_ylabel("AUC (95% CI)", fontsize=fs_axis)
    axD.set_title("Broader temporal discrimination across cohorts", fontsize=fs_panel, fontweight="bold")
    axD.tick_params(labelsize=fs_tick)
    axD.legend(frameon=False, fontsize=fs_note, loc="upper right", bbox_to_anchor=(0.90, 1.0))

    for ax, label in zip([axA, axB, axC, axD], ["A", "B", "C", "D"]):
        ax.text(-0.12, 1.09, label, transform=ax.transAxes, fontsize=max(13 * scale, 10), fontweight="bold")
    fig.suptitle("Generalizability of XGBoost-guided risk stratification", fontsize=fs_title, fontweight="bold", y=0.99)
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.08, right=0.98)
    p1 = safe_savefig(fig, out_png, dpi=350, bbox_inches="tight")
    p2 = safe_savefig(fig, out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return p1, p2


def main():
    root = Path(r"e:\ZhouFX")
    features, ens_base = load_top_features(root, n=120)
    tcga_df = load_tcga(root)
    gse_df = load_gse31210(root, features, ens_base)
    tcga_scores, gse_scores = train_predict_scores(tcga_df, gse_df, features)
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    _, out_main = draw_figure7(
        tcga_df,
        tcga_scores,
        gse_df,
        gse_scores,
        out_dir / "Figure7_rebuilt.png",
        out_dir / "Figure7_rebuilt.pdf",
        figsize=(13.2, 9.8),
        scale=1.0,
    )
    _, out_2col = draw_figure7(
        tcga_df,
        tcga_scores,
        gse_df,
        gse_scores,
        out_dir / "Figure7_rebuilt_2col.png",
        out_dir / "Figure7_rebuilt_2col.pdf",
        figsize=(7.2, 5.9),
        scale=0.62,
    )
    final_main = root / "投稿文件" / "Figure7.pdf"
    shutil.copy2(out_main, final_main)
    print(out_main)
    print(out_2col)
    print(final_main)


if __name__ == "__main__":
    main()
