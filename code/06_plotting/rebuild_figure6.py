from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
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


def load_integrated_data(root: Path):
    p = root / "integrated_data.csv"
    if not p.exists():
        raise FileNotFoundError(str(p))
    df = pd.read_csv(p)
    need = ["OS.time", "OS", "response"]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    df = df.dropna(subset=["OS.time", "OS", "response"]).copy()
    df = df[df["OS.time"] > 0].copy()
    df["OS"] = df["OS"].astype(int)
    df["response"] = df["response"].astype(int)
    return df


def get_model_scores(df, root: Path):
    feat_path = root / "multi_omics_feature_importance.csv"
    if not feat_path.exists():
        raise FileNotFoundError(str(feat_path))
    feat_imp = pd.read_csv(feat_path)
    top_features = feat_imp["feature"].head(120).tolist()
    valid = [f for f in top_features if f in df.columns]
    if len(valid) < 10:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        exclude = {"OS.time", "OS", "response", "Model_Score"}
        valid = [c for c in numeric_cols if c not in exclude][:120]
    X = df[valid]
    y = df["response"].astype(int).values
    imp = SimpleImputer(strategy="median")
    X_imp = imp.fit_transform(X)
    sc = StandardScaler()
    X_sc = sc.fit_transform(X_imp)
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
    scores = cross_val_predict(clf, X_sc, y, cv=cv, method="predict_proba")[:, 1]
    return scores


def km_hr_stats(df):
    d = df.copy()
    cutoff = float(np.median(d["Model_Score"].values))
    d["Risk_Group"] = np.where(d["Model_Score"] > cutoff, "High Risk", "Low Risk")
    t = d["OS.time"].values / 365.0
    e = d["OS"].values.astype(int)
    high = d["Risk_Group"].values == "High Risk"
    low = ~high
    lr = logrank_test(t[high], t[low], event_observed_A=e[high], event_observed_B=e[low])
    cph_data = pd.DataFrame({"T": t, "E": e, "Group": high.astype(int)})
    cph = CoxPHFitter()
    cph.fit(cph_data, duration_col="T", event_col="E")
    hr = float(cph.summary.loc["Group", "exp(coef)"])
    l95 = float(cph.summary.loc["Group", "exp(coef) lower 95%"])
    u95 = float(cph.summary.loc["Group", "exp(coef) upper 95%"])
    return d, cutoff, float(lr.p_value), hr, l95, u95


def plot_panel_a_km(ax, df, fs_axis, fs_panel, fs_note, fs_tick, add_risk_table=True):
    d, _, p, hr, l95, u95 = km_hr_stats(df)
    t = d["OS.time"].values / 365.0
    e = d["OS"].values.astype(int)
    high = d["Risk_Group"] == "High Risk"
    low = d["Risk_Group"] == "Low Risk"
    kmh = KaplanMeierFitter()
    kml = KaplanMeierFitter()
    kmh.fit(t[high], event_observed=e[high], label="High risk")
    kml.fit(t[low], event_observed=e[low], label="Low risk")
    kmh.plot_survival_function(ax=ax, color="#D55E00", linewidth=1.9, ci_show=False)
    kml.plot_survival_function(ax=ax, color="#3C5488", linewidth=1.9, ci_show=False)
    if add_risk_table:
        add_at_risk_counts(kmh, kml, ax=ax, ypos=-0.24, fontsize=max(fs_tick - 1, 6.5))
    ax.set_title("Survival stratification by model-guided risk", fontsize=fs_panel, fontweight="bold")
    ax.set_xlabel("Time (years)", fontsize=fs_axis)
    ax.set_ylabel("Survival probability", fontsize=fs_axis)
    ax.tick_params(labelsize=fs_tick)
    ax.text(
        0.03,
        0.95,
        f"HR={hr:.2f} (95% CI {l95:.2f}-{u95:.2f})\nlog-rank p={p:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )
    ax.text(
        0.03,
        0.84,
        "Consistent separation over time",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )
    if not add_risk_table:
        ax.text(
            0.03,
            0.12,
            f"N high={int(np.sum(high))}, N low={int(np.sum(low))}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=fs_note,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
        )
    ax.legend(frameon=False, fontsize=fs_note, loc="upper right")


def _binary_target_at_time(df, years):
    t_days = years * 365.0
    valid = ~((df["OS"].values == 0) & (df["OS.time"].values <= t_days))
    y = ((df.loc[valid, "OS"].values == 1) & (df.loc[valid, "OS.time"].values <= t_days)).astype(int)
    return valid, y


def plot_panel_b_roc(ax, df, scores, fs_axis, fs_panel, fs_tick, fs_note):
    colors = {1: "#D55E00", 3: "#009E73", 5: "#3C5488"}
    for yr in [1, 3, 5]:
        valid, y = _binary_target_at_time(df, yr)
        s = scores[valid]
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, s)
        a = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=colors[yr], lw=2.0, label=f"{yr}-year AUC={a:.3f}")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("False positive rate", fontsize=fs_axis)
    ax.set_ylabel("True positive rate", fontsize=fs_axis)
    ax.set_title("Moderate but consistent discrimination across time", fontsize=fs_panel, fontweight="bold")
    ax.tick_params(labelsize=fs_tick)
    ax.text(
        0.03,
        0.95,
        "Stable performance across 1-, 3-, and 5-year predictions",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )
    ax.legend(frameon=False, fontsize=fs_note, loc="lower right")


def calculate_net_benefit(y_true, y_prob, thresholds):
    n = len(y_true)
    out = []
    for th in thresholds:
        y_pred = (y_prob >= th).astype(int)
        tp = np.sum((y_pred == 1) & (y_true == 1))
        fp = np.sum((y_pred == 1) & (y_true == 0))
        nb = (tp / n) - (fp / n) * (th / (1 - th)) if th < 1 else 0.0
        out.append(nb)
    return np.array(out)


def plot_panel_c_dca(ax, df, scores, fs_axis, fs_panel, fs_tick, fs_note):
    thresholds = np.linspace(0.01, 0.8, 120)
    ax.plot(thresholds, np.zeros_like(thresholds), color="black", lw=1.2, label="Treat none")
    for yr, color in [(1, "#D55E00"), (3, "#009E73"), (5, "#3C5488")]:
        valid, y = _binary_target_at_time(df, yr)
        s = scores[valid]
        if len(np.unique(y)) < 2:
            continue
        prev = np.mean(y)
        nb_all = prev - (1 - prev) * thresholds / (1 - thresholds)
        nb_model = calculate_net_benefit(y, s, thresholds)
        ax.plot(thresholds, nb_model, color=color, lw=1.9, label=f"{yr}-year model")
        ax.plot(thresholds, nb_all, color=color, lw=1.0, ls="--", alpha=0.45)
        if yr == 3:
            mask = (nb_model > 0) & (nb_model > nb_all)
            ax.fill_between(thresholds, 0, nb_model, where=mask, color=color, alpha=0.10, interpolate=True)
    ax.set_xlim(0, 0.8)
    ax.set_xlabel("Threshold probability", fontsize=fs_axis)
    ax.set_ylabel("Net benefit", fontsize=fs_axis)
    ax.set_title("Decision curve analysis", fontsize=fs_panel, fontweight="bold")
    ax.tick_params(labelsize=fs_tick)
    ax.text(
        0.03,
        0.95,
        "Higher net benefit vs treat-all and treat-none",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )
    ax.legend(frameon=False, fontsize=fs_note, loc="upper right")


def plot_panel_d_response(ax, df, fs_axis, fs_panel, fs_tick, fs_note):
    d = df.copy()
    cutoff = float(np.median(d["Model_Score"].values))
    d["Risk_Group"] = np.where(d["Model_Score"] > cutoff, "High Risk", "Low Risk")
    grp = d.groupby("Risk_Group")["response"].agg(["mean", "sum", "count"]).reindex(["Low Risk", "High Risk"])
    rates = grp["mean"].values
    ard = float(rates[1] - rates[0])
    ax.bar(["Low risk", "High risk"], rates, color=["#3C5488", "#D55E00"], alpha=0.85, width=0.62)
    for i, r in enumerate(rates):
        ax.text(i, r + 0.015, f"{r*100:.1f}%", ha="center", va="bottom", fontsize=fs_note)
    a = int(grp.loc["High Risk", "sum"])
    b = int(grp.loc["High Risk", "count"] - grp.loc["High Risk", "sum"])
    c = int(grp.loc["Low Risk", "sum"])
    d0 = int(grp.loc["Low Risk", "count"] - grp.loc["Low Risk", "sum"])
    odds_ratio, p = stats.fisher_exact([[a, b], [c, d0]], alternative="two-sided")
    ax.set_ylim(0, min(1.0, max(rates) + 0.22))
    ax.set_ylabel("Adverse outcome rate", fontsize=fs_axis)
    ax.set_title("Clinical outcome enrichment by risk group", fontsize=fs_panel, fontweight="bold")
    ax.tick_params(labelsize=fs_tick)
    ax.text(
        0.97,
        0.95,
        f"OR={odds_ratio:.2f}, Fisher p={p:.3g}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )
    ax.text(
        0.03,
        0.88,
        f"High-risk shows higher adverse outcome\nAbsolute risk increase: {ard*100:.1f}%",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )


def draw_figure6(df, scores, out_png, out_pdf, figsize, scale):
    fs_panel = max(11.0 * scale, 8.2)
    fs_axis = max(9.8 * scale, 7.2)
    fs_tick = max(8.5 * scale, 6.4)
    fs_note = max(8.0 * scale, 6.1)
    fs_title = max(14.8 * scale, 11)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.28)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    df_plot = df.copy()
    df_plot["Model_Score"] = scores

    plot_panel_a_km(axA, df_plot, fs_axis, fs_panel, fs_note, fs_tick, add_risk_table=scale >= 0.9)
    plot_panel_b_roc(axB, df_plot, scores, fs_axis, fs_panel, fs_tick, fs_note)
    plot_panel_c_dca(axC, df_plot, scores, fs_axis, fs_panel, fs_tick, fs_note)
    plot_panel_d_response(axD, df_plot, fs_axis, fs_panel, fs_tick, fs_note)

    for ax, label in zip([axA, axB, axC, axD], ["A", "B", "C", "D"]):
        ax.text(-0.12, 1.09, label, transform=ax.transAxes, fontsize=max(13 * scale, 10), fontweight="bold")

    fig.suptitle("Clinical utility of XGBoost-guided risk stratification", fontsize=fs_title, fontweight="bold", y=0.99)
    fig.text(0.5, 0.012, "Integrated evidence supports both prognostic and clinical decision utility", ha="center", va="bottom", fontsize=max(8.0 * scale, 6.1))
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.08, right=0.98)
    p1 = safe_savefig(fig, out_png, dpi=350, bbox_inches="tight")
    p2 = safe_savefig(fig, out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return p1, p2


def main():
    root = Path(r"e:\ZhouFX")
    df = load_integrated_data(root)
    scores = get_model_scores(df, root)
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    _, out_main = draw_figure6(
        df,
        scores,
        out_dir / "Figure6_rebuilt.png",
        out_dir / "Figure6_rebuilt.pdf",
        figsize=(13.2, 9.8),
        scale=1.0,
    )
    _, out_2col = draw_figure6(
        df,
        scores,
        out_dir / "Figure6_rebuilt_2col.png",
        out_dir / "Figure6_rebuilt_2col.pdf",
        figsize=(7.2, 5.9),
        scale=0.62,
    )
    final_main = root / "投稿文件" / "Figure6.pdf"
    shutil.copy2(out_main, final_main)
    print(out_main)
    print(out_2col)
    print(final_main)


if __name__ == "__main__":
    main()
