from pathlib import Path
import shutil
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from scipy import sparse
from scipy.stats import pearsonr
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
        if coords.shape[1] >= 2:
            return coords[:, 0].astype(float), coords[:, 1].astype(float)
    if "array_col" in adata.obs.columns and "array_row" in adata.obs.columns:
        return adata.obs["array_col"].astype(float).values, adata.obs["array_row"].astype(float).values
    raise ValueError("No spatial coordinates in adata.obsm['spatial'] or adata.obs[array_col/array_row]")


def build_mechanism_table(adata):
    x, y = get_spatial_xy(adata)
    c1qa = get_gene_vector(adata, "C1QA")
    c1qb = get_gene_vector(adata, "C1QB")
    c1qc = get_gene_vector(adata, "C1QC")
    spp1 = get_gene_vector(adata, "SPP1")
    cd8a = get_gene_vector(adata, "CD8A")
    cd8b = get_gene_vector(adata, "CD8B")
    c1q = np.mean(np.vstack([c1qa, c1qb, c1qc]), axis=0)
    cd8 = np.mean(np.vstack([cd8a, cd8b]), axis=0) if np.any(cd8b) else cd8a
    df = pd.DataFrame({"x": x, "y": y, "c1q": c1q, "cd8": cd8, "spp1": spp1})
    return df


def load_cosmx_flat_table(root: Path, sample_size=None, seed=123):
    base = root / "Data" / "All+SMI+Flat+data" / "Lung12" / "Lung12-Flat_files_and_images"
    expr = pd.read_csv(base / "Lung12_exprMat_file.csv")
    meta = pd.read_csv(base / "Lung12_metadata_file.csv")
    gene_cols = [c for c in expr.columns if c not in {"fov", "cell_ID"}]
    expr = expr[["cell_ID"] + gene_cols].copy()
    expr["cell_ID"] = expr["cell_ID"].astype(str)
    meta["cell_ID"] = meta["cell_ID"].astype(str)
    merged = meta.merge(expr, on="cell_ID", how="inner")
    c1qa = merged["C1QA"].values if "C1QA" in merged.columns else np.zeros(len(merged))
    c1qb = merged["C1QB"].values if "C1QB" in merged.columns else np.zeros(len(merged))
    c1qc = merged["C1QC"].values if "C1QC" in merged.columns else np.zeros(len(merged))
    cd8a = merged["CD8A"].values if "CD8A" in merged.columns else np.zeros(len(merged))
    cd8b = merged["CD8B"].values if "CD8B" in merged.columns else np.zeros(len(merged))
    spp1 = merged["SPP1"].values if "SPP1" in merged.columns else np.zeros(len(merged))
    c1q = np.mean(np.vstack([c1qa, c1qb, c1qc]), axis=0)
    cd8 = np.mean(np.vstack([cd8a, cd8b]), axis=0) if np.any(cd8b) else cd8a
    q_c1q = np.quantile(c1q, 0.9)
    q_cd8 = np.quantile(cd8, 0.9)
    q_spp1 = np.quantile(spp1, 0.9)
    cell_type = np.array(["Other"] * len(merged), dtype=object)
    cell_type[c1q > q_c1q] = "C1Q+_Macrophage"
    cell_type[cd8 > q_cd8] = "CD8+_T_cell"
    cell_type[spp1 > q_spp1] = "SPP1+_Cell"
    out = pd.DataFrame(
        {
            "cell_ID": merged["cell_ID"].values,
            "x": merged["CenterX_global_px"].astype(float).values,
            "y": merged["CenterY_global_px"].astype(float).values,
            "c1q": c1q.astype(float),
            "cd8": cd8.astype(float),
            "spp1": spp1.astype(float),
            "cell_type_def": cell_type,
        }
    )
    if sample_size is not None and sample_size > 0 and len(out) > sample_size:
        out = out.sample(n=sample_size, random_state=seed).reset_index(drop=True)
    return out


def run_neighborhood_permutation(df, n_permutations=1000):
    coords = df[["x", "y"]].values
    labels = df["cell_type_def"].astype(str).values
    n_neighbors = min(7, len(coords))
    nn = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
    _, idx = nn.kneighbors(coords)
    edges = set()
    for i in range(len(coords)):
        for j in idx[i, 1:]:
            a, b = (i, int(j)) if i < int(j) else (int(j), i)
            if a != b:
                edges.add((a, b))
    edges = np.array(sorted(edges), dtype=int)
    edge_total = float(len(edges))
    targets = ["C1Q+_Macrophage--CD8+_T_cell", "C1Q+_Macrophage--SPP1+_Cell", "CD8+_T_cell--SPP1+_Cell"]
    first = labels[edges[:, 0]]
    second = labels[edges[:, 1]]
    obs = {}
    for t in targets:
        a, b = t.split("--")
        obs[t] = float((((first == a) & (second == b)) | ((first == b) & (second == a))).sum()) / edge_total
    rng = np.random.default_rng(123)
    sims = {t: np.zeros(n_permutations, dtype=float) for t in targets}
    for i in range(n_permutations):
        perm = rng.permutation(labels)
        p1 = perm[edges[:, 0]]
        p2 = perm[edges[:, 1]]
        for t in targets:
            a, b = t.split("--")
            sims[t][i] = float((((p1 == a) & (p2 == b)) | ((p1 == b) & (p2 == a))).sum()) / edge_total
    out = []
    for t in targets:
        s = sims[t]
        mu = float(np.mean(s))
        sd = float(np.std(s, ddof=1))
        z = (obs[t] - mu) / sd if sd > 0 else 0.0
        p = (np.sum(np.abs(s - mu) >= abs(obs[t] - mu)) + 1) / (n_permutations + 1)
        out.append({"unified_int": t, "observed": obs[t], "expected": mu, "z_score": z, "p_value": float(p), "n_permutations": int(n_permutations)})
    return pd.DataFrame(out)


def collect_core_metrics(df, ndf):
    hi_c1q = df["cell_type_def"].values == "C1Q+_Macrophage"
    hi_cd8 = df["cell_type_def"].values == "CD8+_T_cell"
    hi_spp1 = df["cell_type_def"].values == "SPP1+_Cell"
    median_dist = np.nan
    if hi_c1q.sum() > 2 and hi_cd8.sum() > 2:
        nbrs = NearestNeighbors(n_neighbors=1).fit(df.loc[hi_cd8, ["x", "y"]].values)
        dist, _ = nbrs.kneighbors(df.loc[hi_c1q, ["x", "y"]].values)
        median_dist = float(np.median(dist))
    overlap = int(np.sum(hi_c1q & hi_spp1))
    union = int(np.sum(hi_c1q | hi_spp1))
    jaccard = overlap / union if union > 0 else np.nan
    r, p = pearsonr(df["c1q"].values, df["spp1"].values)
    metrics = {
        "median_nearest_distance_c1q_to_cd8": median_dist,
        "r_c1q_spp1": float(r),
        "p_c1q_spp1": float(p),
        "jaccard_c1q_spp1": float(jaccard),
    }
    for _, row in ndf.iterrows():
        key = row["unified_int"].replace("+", "plus").replace("-", "_").replace("__", "_")
        metrics[f"z_{key}"] = float(row["z_score"])
        metrics[f"p_{key}"] = float(row["p_value"])
    return metrics


def _sample_for_display(df, type_col, type_name, max_n, seed):
    sub = df[df[type_col] == type_name]
    if len(sub) > max_n:
        return sub.sample(n=max_n, random_state=seed)
    return sub


def run_subsampling_robustness(root: Path, sample_size=12000, seeds=(1, 2, 3), n_permutations=1000):
    rows = []
    for s in seeds:
        df = load_cosmx_flat_table(root, sample_size=sample_size, seed=s)
        ndf = run_neighborhood_permutation(df, n_permutations=n_permutations)
        rec = {"seed": int(s), "sample_size": int(sample_size), "n_permutations": int(n_permutations)}
        rec.update(collect_core_metrics(df, ndf))
        rows.append(rec)
    out = pd.DataFrame(rows)
    out_path = root / "Figure4_subsampling_robustness.csv"
    out.to_csv(out_path, index=False)
    return out_path


def panel_a_colocalization(ax, df):
    bg = df if len(df) <= 3500 else df.sample(n=3500, random_state=101)
    c1q_df = _sample_for_display(df, "cell_type_def", "C1Q+_Macrophage", 900, 102)
    cd8_df = _sample_for_display(df, "cell_type_def", "CD8+_T_cell", 900, 103)
    ax.scatter(bg["x"], bg["y"], s=0.9, c="#D1D5DB", alpha=0.14, linewidths=0)
    ax.scatter(c1q_df["x"], c1q_df["y"], s=5.2, c="#D55E00", alpha=0.72, label="C1Q+ macrophages")
    ax.scatter(
        cd8_df["x"],
        cd8_df["y"],
        s=6.2,
        facecolors="none",
        edgecolors="#009E73",
        linewidths=0.5,
        alpha=0.8,
        label="CD8+ T cells",
    )
    if len(c1q_df) > 2 and len(cd8_df) > 2:
        nbrs = NearestNeighbors(n_neighbors=1).fit(cd8_df[["x", "y"]].values)
        dist, _ = nbrs.kneighbors(c1q_df[["x", "y"]].values)
        median_dist = float(np.median(dist))
        ax.text(
            0.02,
            0.98,
            f"Median distance = {median_dist:.1f}",
            transform=ax.transAxes,
            fontsize=8.2,
            va="top",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
        )
    ax.set_title("Spatial association of C1Q+ macrophages and CD8+ T cells", fontsize=10.5, fontweight="bold")
    ax.set_xlabel("Spatial x")
    ax.set_ylabel("Spatial y")


def panel_b_neighborhood(ax, ndf):
    show = ndf.copy().sort_values("z_score", ascending=True)
    label_map = {
        "C1Q+_Macrophage--CD8+_T_cell": "C1Q+ Mac–CD8+ T",
        "C1Q+_Macrophage--SPP1+_Cell": "C1Q+ Mac–SPP1+",
        "CD8+_T_cell--SPP1+_Cell": "CD8+ T–SPP1+",
    }
    labels = [label_map.get(x, x.replace("_", " ")) for x in show["unified_int"].values]
    vals = show["z_score"].values
    colors = ["#D55E00" if "C1Q+_Macrophage--CD8+_T_cell" in s else "#3C5488" for s in show["unified_int"].values]
    ax.barh(labels, vals, color=colors, alpha=0.9)
    max_abs = max(0.5, float(np.max(np.abs(vals))) if len(vals) else 0.5)
    text_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for i, (_, r) in enumerate(show.iterrows()):
        ax.text(0.88, i, f"p={r['p_value']:.2f}", transform=text_transform, va="center", ha="left", fontsize=8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlim(-max_abs * 1.1, max_abs * 1.95)
    ax.set_xlabel("Permutation z-score (1000 permutations)")
    ax.set_title("Neighborhood enrichment analysis", fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="y", labelsize=8.2, pad=2)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", labelright=True, labelleft=False)


def panel_c_correlation(ax, adata, df):
    markers = ["CD163", "CD68", "CD14", "CD8A", "CTLA4", "GZMB", "SPP1"]
    out = []
    for g in markers:
        if g in adata.var_names:
            v = get_gene_vector(adata, g)
            r, p = pearsonr(df["c1q"].values, v)
            out.append((g, float(r), float(p)))
    cdf = pd.DataFrame(out, columns=["gene", "r", "p"])
    if cdf.empty:
        ax.text(0.5, 0.5, "No marker available", ha="center", va="center")
        ax.set_axis_off()
        return
    cdf = cdf.sort_values("r", ascending=True)
    colors = ["#0072B2" if r < 0 else "#D55E00" for r in cdf["r"]]
    ax.barh(cdf["gene"], cdf["r"], color=colors, alpha=0.9)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlim(min(-0.4, cdf["r"].min() - 0.08), max(0.4, cdf["r"].max() + 0.08))
    ax.set_xlabel("Pearson r with C1Q score")
    ax.set_title("C1Q–immune correlation", fontsize=10.5, fontweight="bold")


def panel_d_antagonism(ax, df):
    hi_c1q = df["cell_type_def"].values == "C1Q+_Macrophage"
    hi_spp1 = df["cell_type_def"].values == "SPP1+_Cell"
    overlap = np.sum(hi_c1q & hi_spp1)
    union = np.sum(hi_c1q | hi_spp1)
    jaccard = overlap / union if union > 0 else 0
    r, p = pearsonr(df["c1q"].values, df["spp1"].values)
    bg = df if len(df) <= 3500 else df.sample(n=3500, random_state=201)
    c1q_df = _sample_for_display(df, "cell_type_def", "C1Q+_Macrophage", 900, 202)
    spp1_df = _sample_for_display(df, "cell_type_def", "SPP1+_Cell", 900, 203)
    ax.scatter(bg["x"], bg["y"], s=0.9, c="#E5E7EB", alpha=0.14, linewidths=0)
    ax.scatter(c1q_df["x"], c1q_df["y"], s=5.2, c="#D55E00", alpha=0.72, label="C1Q+ macrophages")
    ax.scatter(
        spp1_df["x"],
        spp1_df["y"],
        s=6.2,
        facecolors="none",
        edgecolors="#0072B2",
        linewidths=0.5,
        alpha=0.85,
        label="SPP1+ cells",
    )
    ax.text(
        0.02,
        0.98,
        f"r(C1Q, SPP1) = {r:.3f}\nJaccard = {jaccard:.2f}",
        transform=ax.transAxes,
        fontsize=8.2,
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, boxstyle="round,pad=0.2"),
    )
    ax.set_title("Spatial antagonism between C1Q and SPP1 niches", fontsize=10.5, fontweight="bold")
    ax.set_xlabel("Spatial x")
    ax.set_ylabel("Spatial y")


def draw_figure4_mechanism(adata, cosmx_df, ndf, out_png, out_pdf, figsize, scale):
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    panel_a_colocalization(axA, cosmx_df)
    panel_b_neighborhood(axB, ndf)
    panel_c_correlation(axC, adata, build_mechanism_table(adata))
    panel_d_antagonism(axD, cosmx_df)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#D55E00", markeredgecolor="#D55E00", markersize=5, label="C1Q+ macrophages"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#009E73", markersize=5, label="CD8+ T cells"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#0072B2", markersize=5, label="SPP1+ cells"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=3, frameon=False, fontsize=max(8.0 * scale, 7.0))

    fs = max(15 * scale, 11)
    fig.suptitle("Spatial Association Pattern", fontsize=fs, fontweight="bold", y=0.99)
    fig.subplots_adjust(top=0.88, bottom=0.07, left=0.06, right=0.98)
    p1 = safe_savefig(fig, out_png, dpi=350, bbox_inches="tight")
    p2 = safe_savefig(fig, out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return p1, p2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-size", type=int, default=12000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--robustness", action="store_true")
    args = parser.parse_args()

    root = Path(r"e:\ZhouFX")
    data_path = pick_existing_path(
        [
            root / "GigaTIME_Paper_Submission" / "analysis_scripts" / "processed_spatial_data.h5ad",
            root / "processed_spatial_data.h5ad",
        ]
    )
    adata = sc.read_h5ad(str(data_path))
    sample_size = None if args.full else args.sample_size
    cosmx_df = load_cosmx_flat_table(root, sample_size=sample_size, seed=args.seed)
    ndf = run_neighborhood_permutation(cosmx_df, n_permutations=1000)

    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    _, out_main = draw_figure4_mechanism(
        adata,
        cosmx_df,
        ndf,
        out_dir / "Figure4_rebuilt.png",
        out_dir / "Figure4_rebuilt.pdf",
        figsize=(13.2, 9.8),
        scale=1.0,
    )
    _, out_2col = draw_figure4_mechanism(
        adata,
        cosmx_df,
        ndf,
        out_dir / "Figure4_rebuilt_2col.png",
        out_dir / "Figure4_rebuilt_2col.pdf",
        figsize=(7.2, 5.8),
        scale=0.62,
    )

    final_main = root / "投稿文件" / "Figure4.pdf"
    shutil.copy2(out_main, final_main)
    ndf.to_csv(root / "neighborhood_enrichment_results.csv", index=False)
    cosmx_df[["cell_ID", "cell_type_def", "c1q", "cd8", "spp1"]].to_csv(root / "cell_type_definition_mapping.csv", index=False)

    print(final_main)
    print(out_2col)
    print(root / "neighborhood_enrichment_results.csv")
    print(root / "cell_type_definition_mapping.csv")
    if args.robustness:
        rb = run_subsampling_robustness(root, sample_size=args.sample_size, seeds=(1, 2, 3), n_permutations=1000)
        print(rb)


if __name__ == "__main__":
    main()
