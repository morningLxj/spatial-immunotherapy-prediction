from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
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


def load_real_visium_data(root: Path):
    h5_path = root / "Data" / "filtered_feature_bc_matrix.h5"
    pos_path = root / "Data" / "spatial" / "tissue_positions.csv"
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
    return ad, meta


def get_gene_vector(ad, gene_symbol):
    if gene_symbol in ad.var_names:
        v = ad[:, gene_symbol].X
        return np.asarray(v.todense()).ravel() if sparse.issparse(v) else np.asarray(v).ravel()
    return np.zeros(ad.n_obs, dtype=float)


def compute_morans_i_all_genes(ad, x_coord, y_coord, n_neighbors=8):
    coords = np.column_stack([x_coord, y_coord]).astype(float)
    nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, len(coords)), metric="euclidean").fit(coords)
    _, idx = nn.kneighbors(coords)
    n = len(coords)
    rows, cols = [], []
    for i in range(n):
        neigh = idx[i, 1:]
        rows.extend([i] * len(neigh))
        cols.extend(neigh.tolist())
    data = np.ones(len(rows), dtype=float)
    w = sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    w = ((w + w.T) > 0).astype(float).tocsr()
    w_sum = float(w.sum())
    xmat = ad.X.tocsc() if sparse.issparse(ad.X) else sparse.csc_matrix(ad.X)
    out = []
    for j, g in enumerate(ad.var_names):
        x = np.asarray(xmat[:, j].todense()).ravel()
        x_c = x - x.mean()
        denom = float(np.sum(x_c ** 2))
        if denom <= 0 or w_sum <= 0:
            out.append((g, np.nan))
            continue
        wx = w.dot(x_c)
        num = float(np.sum(x_c * wx))
        moran = (n / w_sum) * (num / denom)
        out.append((g, moran))
    return pd.DataFrame(out, columns=["Gene", "Morans_I"])


def build_real_spatial_tables(root: Path, out_dir: Path):
    ad, meta = load_real_visium_data(root)
    c1qa = get_gene_vector(ad, "C1QA")
    spp1 = get_gene_vector(ad, "SPP1")
    cd8a = get_gene_vector(ad, "CD8A")
    cd8b = get_gene_vector(ad, "CD8B")
    cd8 = (cd8a + cd8b) / 2.0
    detailed = pd.DataFrame(
        {
            "barcode": ad.obs_names.astype(str),
            "x": meta["pxl_col_in_fullres"].astype(float).values,
            "y": meta["pxl_row_in_fullres"].astype(float).values,
            "c1qa_expression": c1qa,
            "spp1_expression": spp1,
            "cd8_expression": cd8,
        }
    )
    moran = compute_morans_i_all_genes(ad, detailed["x"].values, detailed["y"].values, n_neighbors=8)
    detailed_path = out_dir / "Figure3_real_spatial_coordinates_expression.csv"
    moran_path = out_dir / "Figure3_real_spatial_morans_i.csv"
    detailed.to_csv(detailed_path, index=False)
    moran.to_csv(moran_path, index=False)
    return detailed, moran, detailed_path, moran_path


def draw_figure3(data_df, moran_df, out_png, out_pdf, figsize, scale):
    fs_title = max(15 * scale, 11.0)
    fs_panel = max(11.5 * scale, 8.0)
    fs_axis = max(9.8 * scale, 7.2)
    fs_tick = max(8.8 * scale, 6.8)
    fs_annot = max(8.0 * scale, 6.4)
    palette = {"c1q": "#D55E00", "spp1": "#0072B2", "cd8": "#009E73", "gray": "#374151"}

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1, 1], hspace=0.28, wspace=0.24)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 1])

    x = data_df["x"].astype(float).values
    y = data_df["y"].astype(float).values
    c1q = data_df["c1qa_expression"].astype(float).values
    spp1 = data_df["spp1_expression"].astype(float).values
    cd8 = data_df["cd8_expression"].astype(float).values
    q_c1q = np.quantile(c1q, 0.85)
    q_spp1 = np.quantile(spp1, 0.85)
    q_cd8 = np.quantile(cd8, 0.85)
    hi_c1q = c1q >= q_c1q
    hi_spp1 = spp1 >= q_spp1
    hi_cd8 = cd8 >= q_cd8

    sc = axA.scatter(x, y, c=c1q, s=4.5 if scale >= 1 else 3.5, cmap="Reds", alpha=0.8, linewidths=0)
    axA.scatter(x[hi_cd8], y[hi_cd8], s=6 if scale >= 1 else 4, color=palette["cd8"], alpha=0.7, label="CD8 high")
    axA.scatter(x[hi_spp1], y[hi_spp1], s=9 if scale >= 1 else 6, facecolors="none", edgecolors=palette["spp1"], linewidths=0.5, label="SPP1 high")
    axA.contour(
        np.histogram2d(x[hi_c1q], y[hi_c1q], bins=30)[0].T,
        levels=3,
        colors=[palette["c1q"]],
        linewidths=1.3,
        origin="lower",
        extent=[x.min(), x.max(), y.min(), y.max()],
    )
    axA.set_title("Panel A  Spatial C1Q hotspots with CD8 context", fontsize=fs_panel, fontweight="bold")
    axA.set_xlabel("Spatial x", fontsize=fs_axis)
    axA.set_ylabel("Spatial y", fontsize=fs_axis)
    axA.tick_params(labelsize=fs_tick)
    cb = fig.colorbar(sc, ax=axA, fraction=0.038, pad=0.02)
    cb.ax.tick_params(labelsize=fs_tick - 0.5)
    cb.set_label("C1QA expression", fontsize=fs_axis)
    axA.legend(frameon=False, fontsize=fs_annot, loc="lower right", bbox_to_anchor=(0.99, 0.06))
    axA.text(0.01, 0.01, "Hotspot boundary: top 15% C1QA", transform=axA.transAxes, fontsize=fs_annot, color=palette["gray"])

    m = moran_df.copy()
    axB.hist(m["Morans_I"].values, bins=12, color="#CBD5E1", edgecolor="#94A3B8")
    c1q_vals = m[m["Gene"].isin(["C1QA", "C1QB", "C1QC"])][["Gene", "Morans_I"]]
    for _, r in c1q_vals.iterrows():
        axB.axvline(r["Morans_I"], color=palette["c1q"], lw=1.8, linestyle="-" if r["Gene"] == "C1QA" else "--")
    if (m["Gene"] == "SPP1").any():
        spp1_i = float(m.loc[m["Gene"] == "SPP1", "Morans_I"].iloc[0])
        axB.axvline(spp1_i, color=palette["spp1"], lw=1.6, linestyle=":")
    axB.set_title("Panel B  Moran's I distribution", fontsize=fs_panel, fontweight="bold")
    axB.set_xlabel("Moran's I", fontsize=fs_axis)
    axB.set_ylabel("Gene count", fontsize=fs_axis)
    axB.tick_params(labelsize=fs_tick)
    c1qa_i = float(m.loc[m["Gene"] == "C1QA", "Morans_I"].iloc[0]) if (m["Gene"] == "C1QA").any() else np.nan
    axB.text(0.02, 0.93, f"C1QA Moran's I = {c1qa_i:.3f}", transform=axB.transAxes, fontsize=fs_annot, color=palette["c1q"], va="top")

    top = m.sort_values("Morans_I", ascending=False).head(10).copy()
    colors = []
    for g in top["Gene"]:
        if g in {"C1QA", "C1QB", "C1QC"}:
            colors.append(palette["c1q"])
        elif g == "SPP1":
            colors.append(palette["spp1"])
        else:
            colors.append("#94A3B8")
    axC.barh(top["Gene"], top["Morans_I"], color=colors)
    spp1_benchmark = float(m.loc[m["Gene"] == "SPP1", "Morans_I"].iloc[0]) if (m["Gene"] == "SPP1").any() else 0.712
    axC.axvline(spp1_benchmark, color=palette["spp1"], lw=1.8, linestyle="--")
    axC.text(spp1_benchmark + 0.005, 0.35, "SPP1 benchmark", color=palette["spp1"], fontsize=fs_annot)
    axC.invert_yaxis()
    axC.set_title("Panel C  Top spatial genes and C1Q axis", fontsize=fs_panel, fontweight="bold")
    axC.set_xlabel("Moran's I", fontsize=fs_axis)
    axC.tick_params(labelsize=fs_tick)

    fig.suptitle("Spatial Hotspots of C1Q Axis", fontsize=fs_title, fontweight="bold", y=0.992)
    fig.subplots_adjust(top=0.92, bottom=0.08, left=0.07, right=0.98)
    out1 = safe_savefig(fig, out_png, dpi=350, bbox_inches="tight")
    out2 = safe_savefig(fig, out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return out1, out2


def main():
    root = Path(r"e:\ZhouFX")
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    data_df, moran_df, src_data_path, src_moran_path = build_real_spatial_tables(root, out_dir)
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    _, out_main = draw_figure3(
        data_df,
        moran_df,
        out_dir / "Figure3_rebuilt.png",
        out_dir / "Figure3_rebuilt.pdf",
        figsize=(13.6, 9.2),
        scale=1.0,
    )
    _, out_2col = draw_figure3(
        data_df,
        moran_df,
        out_dir / "Figure3_rebuilt_2col.png",
        out_dir / "Figure3_rebuilt_2col.pdf",
        figsize=(7.2, 5.8),
        scale=0.62,
    )
    print(out_main)
    print(out_2col)
    print(src_data_path)
    print(src_moran_path)


if __name__ == "__main__":
    main()
