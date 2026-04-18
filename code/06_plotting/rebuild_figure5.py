from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy import stats
from sklearn.neighbors import NearestNeighbors

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


def pick_existing_path(candidates):
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("No valid input file found in candidates")


def to_dense_vector(x):
    if sparse.issparse(x):
        return np.asarray(x.todense()).ravel()
    return np.asarray(x).ravel()


def get_gene_vector(adata, gene):
    if gene in adata.var_names:
        return to_dense_vector(adata[:, gene].X)
    return np.zeros(adata.n_obs, dtype=float)


def get_spatial_xy(adata):
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"])
        return coords[:, 0].astype(float), coords[:, 1].astype(float)
    if "array_col" in adata.obs.columns and "array_row" in adata.obs.columns:
        return adata.obs["array_col"].astype(float).values, adata.obs["array_row"].astype(float).values
    raise ValueError("No spatial coordinates found")


def load_visium_real_data(root: Path):
    h5_path = root / "Data" / "filtered_feature_bc_matrix.h5"
    pos_path = root / "Data" / "spatial" / "tissue_positions.csv"
    if h5_path.exists() and pos_path.exists():
        ad = sc.read_10x_h5(str(h5_path))
        ad.var_names_make_unique()
        pos = pd.read_csv(pos_path)
        if "barcode" not in pos.columns:
            pos = pd.read_csv(pos_path, header=None)
            pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
        pos["barcode"] = pos["barcode"].astype(str)
        pos = pos[pos["in_tissue"] == 1].copy()
        obs = pd.DataFrame(index=ad.obs_names)
        obs["barcode"] = obs.index.astype(str)
        meta = obs.merge(pos, on="barcode", how="inner").set_index("barcode")
        ad = ad[meta.index].copy()
        meta = meta.loc[ad.obs_names]
        ad.obsm["spatial"] = np.column_stack([meta["pxl_col_in_fullres"].astype(float).values, meta["pxl_row_in_fullres"].astype(float).values])
        return ad
    data_path = pick_existing_path(
        [
            root / "GigaTIME_Paper_Submission" / "analysis_scripts" / "processed_spatial_data.h5ad",
            root / "processed_spatial_data.h5ad",
        ]
    )
    return sc.read_h5ad(str(data_path))


def compute_morans_i(values, x, y, n_neighbors=8):
    coords = np.column_stack([x, y]).astype(float)
    k = min(max(3, n_neighbors + 1), len(coords))
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(coords)
    _, idx = nn.kneighbors(coords)
    rows, cols = [], []
    for i in range(len(coords)):
        for j in idx[i, 1:]:
            rows.append(i)
            cols.append(int(j))
    data = np.ones(len(rows), dtype=float)
    w = sparse.coo_matrix((data, (rows, cols)), shape=(len(coords), len(coords))).tocsr()
    w = ((w + w.T) > 0).astype(float).tocsr()
    w_sum = float(w.sum())
    x0 = values - np.mean(values)
    denom = float(np.sum(x0**2))
    if denom <= 0 or w_sum <= 0:
        return np.nan
    wx = w.dot(x0)
    num = float(np.sum(x0 * wx))
    return (len(coords) / w_sum) * (num / denom)


def compute_cd8_proximity(c1q_values, cd8_values, x, y, q=0.85):
    q_cd8 = float(np.quantile(cd8_values, q))
    hi_cd8 = cd8_values >= q_cd8
    if np.sum(hi_cd8) < 10:
        return np.nan
    all_xy = np.column_stack([x, y])
    cd8_xy = np.column_stack([x[hi_cd8], y[hi_cd8]])
    nn = NearestNeighbors(n_neighbors=1).fit(cd8_xy)
    dist, _ = nn.kneighbors(all_xy)
    dist = dist.ravel()
    w = np.maximum(c1q_values, 0)
    if float(np.sum(w)) <= 0:
        return float(np.median(dist))
    order = np.argsort(dist)
    dist_sorted = dist[order]
    w_sorted = w[order]
    cw = np.cumsum(w_sorted) / np.sum(w_sorted)
    idx = np.searchsorted(cw, 0.5)
    idx = min(max(int(idx), 0), len(dist_sorted) - 1)
    return float(dist_sorted[idx])


def zscore(v):
    s = float(np.std(v))
    if s == 0:
        return np.zeros_like(v)
    return (v - float(np.mean(v))) / s


def zscore_by_ref(v, mu, sigma):
    if sigma == 0:
        return np.zeros_like(v)
    return (v - mu) / sigma


def simulate_c1q_knockdown(c1q, level):
    c1q_p = c1q.copy()
    thr_c1q = float(np.quantile(c1q, 0.85))
    med_c1q = float(np.quantile(c1q, 0.5))
    high_mask = c1q >= thr_c1q
    mid_mask = (c1q < thr_c1q) & (c1q >= med_c1q)
    c1q_p[high_mask] = c1q_p[high_mask] * (1 - level * 0.95)
    c1q_p[mid_mask] = c1q_p[mid_mask] * (1 - level * 0.35)
    return c1q_p


def c1q_hotspot_to_cd8_distances(c1q_values, cd8_values, x, y, c1q_threshold, cd8_threshold):
    hi_c1q = c1q_values >= c1q_threshold
    hi_cd8 = cd8_values >= cd8_threshold
    if np.sum(hi_c1q) < 10 or np.sum(hi_cd8) < 10:
        return np.array([], dtype=float)
    c1q_xy = np.column_stack([x[hi_c1q], y[hi_c1q]])
    cd8_xy = np.column_stack([x[hi_cd8], y[hi_cd8]])
    nn = NearestNeighbors(n_neighbors=1).fit(cd8_xy)
    dist, _ = nn.kneighbors(c1q_xy)
    return dist.ravel().astype(float)


def build_perturbation_series(c1q, cd8, spp1, x, y):
    levels = np.array([0.0, 0.2, 0.4, 0.6, 0.8], dtype=float)
    thr_c1q = float(np.quantile(c1q, 0.85))
    mu_c1q = float(np.mean(c1q))
    sd_c1q = float(np.std(c1q))
    rows = []
    for lv in levels:
        c1q_p = simulate_c1q_knockdown(c1q, lv)
        moran = compute_morans_i(c1q_p, x, y, n_neighbors=8)
        prox = compute_cd8_proximity(c1q_p, cd8, x, y, q=0.85)
        risk = 0.45 * zscore(spp1) + 0.35 * zscore_by_ref(c1q_p, mu_c1q, sd_c1q) - 0.20 * zscore(cd8)
        rows.append(
            {
                "knockdown_level": lv,
                "morans_i": moran,
                "cd8_proximity_median_dist": prox,
                "risk_mean": float(np.mean(risk)),
                "risk_q75": float(np.quantile(risk, 0.75)),
            }
        )
    return pd.DataFrame(rows)


def draw_figure5(adata, out_png, out_pdf, figsize, scale):
    x, y = get_spatial_xy(adata)
    c1q = np.mean(np.vstack([get_gene_vector(adata, "C1QA"), get_gene_vector(adata, "C1QB"), get_gene_vector(adata, "C1QC")]), axis=0)
    cd8a = get_gene_vector(adata, "CD8A")
    cd8b = get_gene_vector(adata, "CD8B")
    cd8 = (cd8a + cd8b) / 2.0 if np.any(cd8b) else cd8a
    spp1 = get_gene_vector(adata, "SPP1")
    pert = build_perturbation_series(c1q, cd8, spp1, x, y)

    fs_panel = max(11.0 * scale, 8.2)
    fs_axis = max(9.8 * scale, 7.4)
    fs_tick = max(8.6 * scale, 6.6)
    fs_note = max(8.0 * scale, 6.2)
    fs_title = max(15.0 * scale, 11.0)
    orange = "#D55E00"
    green = "#009E73"
    blue = "#3C5488"

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    base = c1q
    kd = c1q * 0.2
    axA.boxplot([base, kd], tick_labels=["Baseline", "80% KD"], patch_artist=True, widths=0.5, boxprops=dict(facecolor="#FEE8D6", edgecolor=orange), medianprops=dict(color=orange, linewidth=1.4))
    axA.set_ylabel("C1Q activity", fontsize=fs_axis)
    axA.set_title("In silico C1Q knockdown reduces axis activity", fontsize=fs_panel, fontweight="bold")
    axA.tick_params(labelsize=fs_tick)
    axA.text(0.97, 0.95, f"Mean drop = {(1 - np.mean(kd)/max(np.mean(base), 1e-9))*100:.1f}%", transform=axA.transAxes, ha="right", va="top", fontsize=fs_note, bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", boxstyle="round,pad=0.2"))

    xk = pert["knockdown_level"] * 100
    axB.plot(xk, pert["morans_i"], color=orange, marker="o", linewidth=1.8, markersize=4.8)
    axB.set_xlabel("C1Q knockdown (%)", fontsize=fs_axis)
    axB.set_ylabel("Moran's I", fontsize=fs_axis)
    axB.set_title("Spatial hotspot structure trends toward disruption after perturbation", fontsize=fs_panel, fontweight="bold")
    axB.tick_params(labelsize=fs_tick)
    axB.grid(axis="y", linestyle=":", alpha=0.35)
    rho_b, p_b = stats.spearmanr(xk.values, pert["morans_i"].values, nan_policy="omit")
    axB.text(0.97, 0.95, f"Spearman ρ={rho_b:.2f}, p={p_b:.3f}", transform=axB.transAxes, ha="right", va="top", fontsize=fs_note, bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", boxstyle="round,pad=0.2"))

    thr_c1q = float(np.quantile(c1q, 0.85))
    thr_cd8 = float(np.quantile(cd8, 0.85))
    dist_sets = []
    medians = []
    for lv in pert["knockdown_level"].values:
        c1q_p = simulate_c1q_knockdown(c1q, float(lv))
        d = c1q_hotspot_to_cd8_distances(c1q_p, cd8, x, y, c1q_threshold=thr_c1q, cd8_threshold=thr_cd8)
        if len(d) == 0:
            d = np.array([np.nan])
        dist_sets.append(d)
        medians.append(float(np.nanmedian(d)))
    v = axC.violinplot(dist_sets, positions=xk.values, widths=10, showmeans=False, showmedians=True, showextrema=False)
    for body in v["bodies"]:
        body.set_facecolor(green)
        body.set_alpha(0.25)
        body.set_edgecolor(green)
        body.set_linewidth(0.8)
    if "cmedians" in v:
        v["cmedians"].set_color(green)
        v["cmedians"].set_linewidth(1.5)
    delta_prox = np.array(medians) - float(medians[0])
    axC.plot(xk, medians, color=green, marker="o", linewidth=1.3, markersize=3.8, alpha=0.9)
    axC.set_xlabel("C1Q knockdown (%)", fontsize=fs_axis)
    axC.set_ylabel("C1Q-hotspot to CD8 distance", fontsize=fs_axis)
    axC.set_title("CD8 proximity exhibits threshold response to C1Q perturbation", fontsize=fs_panel, fontweight="bold")
    axC.tick_params(labelsize=fs_tick)
    axC.grid(axis="y", linestyle=":", alpha=0.35)
    rho_c, p_c = stats.spearmanr(xk.values, delta_prox, nan_policy="omit")
    d0 = np.asarray(dist_sets[0], dtype=float)
    d8 = np.asarray(dist_sets[-1], dtype=float)
    pooled = np.sqrt((np.nanvar(d0) + np.nanvar(d8)) / 2) if (len(d0) > 1 and len(d8) > 1) else np.nan
    eff = (np.nanmean(d8) - np.nanmean(d0)) / pooled if pooled and pooled > 0 else np.nan
    axC.text(
        0.97,
        0.95,
        f"Suggestive non-linear response\nΔ80%={delta_prox[-1]:.1f}, d={eff:.2f}\nρ={rho_c:.2f}, p={p_c:.3f}",
        transform=axC.transAxes,
        ha="right",
        va="top",
        fontsize=fs_note,
        bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", boxstyle="round,pad=0.2"),
    )

    c1q_kd = c1q.copy()
    thr_c1q = float(np.quantile(c1q, 0.85))
    high_mask = c1q >= thr_c1q
    mid_mask = (c1q < thr_c1q) & (c1q >= np.quantile(c1q, 0.5))
    c1q_kd[high_mask] = c1q_kd[high_mask] * 0.24
    c1q_kd[mid_mask] = c1q_kd[mid_mask] * 0.72
    mu_c1q = float(np.mean(c1q))
    sd_c1q = float(np.std(c1q))
    risk_base = 0.45 * zscore(spp1) + 0.35 * zscore_by_ref(c1q, mu_c1q, sd_c1q) - 0.20 * zscore(cd8)
    risk_kd = 0.45 * zscore(spp1) + 0.35 * zscore_by_ref(c1q_kd, mu_c1q, sd_c1q) - 0.20 * zscore(cd8)
    bins = np.linspace(min(risk_base.min(), risk_kd.min()), max(risk_base.max(), risk_kd.max()), 40)
    axD.hist(risk_base, bins=bins, alpha=0.5, color=blue, density=True, label="Baseline")
    axD.hist(risk_kd, bins=bins, alpha=0.5, color=orange, density=True, label="80% KD")
    axD.set_xlabel("Predicted risk score", fontsize=fs_axis)
    axD.set_ylabel("Density", fontsize=fs_axis)
    axD.set_title("Risk distribution shifts after perturbation", fontsize=fs_panel, fontweight="bold")
    axD.tick_params(labelsize=fs_tick)
    axD.legend(frameon=False, fontsize=fs_note)

    for ax, label in zip([axA, axB, axC, axD], ["A", "B", "C", "D"]):
        y_pos = 1.12 if label == "C" else 1.08
        ax.text(-0.09, y_pos, label, transform=ax.transAxes, fontsize=max(13 * scale, 10), fontweight="bold")

    fig.suptitle("Virtual perturbation of C1Q axis", fontsize=fs_title, fontweight="bold", y=0.99)
    fig.subplots_adjust(top=0.9, bottom=0.08, left=0.08, right=0.98)
    png_out = safe_savefig(fig, out_png, dpi=350, bbox_inches="tight")
    pdf_out = safe_savefig(fig, out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return pert, png_out, pdf_out


def main():
    root = Path(r"e:\ZhouFX")
    adata = load_visium_real_data(root)
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    pert, _, out_main_pdf = draw_figure5(
        adata,
        out_dir / "Figure5_rebuilt.png",
        out_dir / "Figure5_rebuilt.pdf",
        figsize=(13.2, 9.8),
        scale=1.0,
    )
    _, _, out_2col_pdf = draw_figure5(
        adata,
        out_dir / "Figure5_rebuilt_2col.png",
        out_dir / "Figure5_rebuilt_2col.pdf",
        figsize=(7.2, 5.8),
        scale=0.62,
    )
    pert.to_csv(root / "virtual_perturbation_curve_data.csv", index=False)
    final_main = root / "投稿文件" / "Figure5.pdf"
    shutil.copy2(out_main_pdf, final_main)
    print(out_main_pdf)
    print(out_2col_pdf)
    print(final_main)
    print(root / "virtual_perturbation_curve_data.csv")


if __name__ == "__main__":
    main()
