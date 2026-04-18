from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches

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


def draw_evidence_driven_figure1(out_png, out_pdf, figsize, scale):
    compact = scale < 0.8
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 1], hspace=0.18, wspace=0.1)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    for ax in [ax_a, ax_b, ax_c]:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    fs_title = max((22 if not compact else 20) * scale, 13.5)
    fs_panel = max((15 if not compact else 13) * scale, 9.0)
    fs_node = max((12.2 if not compact else 10.8) * scale, 8.2)
    fs_arrow = max((10.1 if not compact else 8.9) * scale, 7.3)
    fs_small = max((9.1 if not compact else 7.8) * scale, 7.0)
    fs_layer = max((11.2 if not compact else 9.8) * scale, 8.0)
    palette = {
        "line": "#1F2937",
        "accent": "#B2182B",
        "c1q": "#D55E00",
        "spp1": "#0072B2",
        "cd8": "#009E73",
        "neutral": "#F8FAFC",
        "node1": "#F4F4F5",
        "node2": "#FFF3E8",
        "node3": "#EDF7ED",
        "node4": "#EAF4FD",
        "node5": "#F7ECFA",
    }
    fig.suptitle("Causal–Spatial–Functional axis of C1Q in NSCLC", fontsize=fs_title, fontweight="bold", y=0.975)
    ax_a.text(0.01, 0.93, "Panel A  Mechanism axis", fontsize=fs_panel, fontweight="bold")
    node_y = 0.42 if not compact else 0.41
    node_h = 0.28 if not compact else 0.21
    nodes = [
        (0.03, node_y, 0.16, node_h, "Genetic\nvariation", palette["node1"]),
        (0.24, node_y, 0.16, node_h, "C1Q\nactivity", palette["node2"]),
        (0.45, node_y, 0.16, node_h, "Spatial\nhotspot", palette["node3"]),
        (0.66, node_y, 0.16, node_h, "CD8\nproximity", palette["node4"]),
        (0.85, node_y, 0.12, node_h, "ICB\nresponse", palette["node5"]),
    ]
    for x, y, w, h, text, color in nodes:
        box = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018", linewidth=1.9, edgecolor=palette["line"], facecolor=color)
        ax_a.add_patch(box)
        ax_a.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs_node, fontweight="bold")
    arrow_y = node_y + node_h * 0.50
    arrow_labels = [
        "MR (P<1e-10)\nSteiger OK",
        "Moran I = 0.66",
        "Neighbor ↑\nMφ ↔ CD8+",
        "Mediation OK",
    ]
    label_tiers = [0.032, 0.014, 0.032, 0.014] if not compact else [0.078, 0.058, 0.078, 0.058]
    for i in range(4):
        x1 = nodes[i][0] + nodes[i][2]
        x2 = nodes[i + 1][0]
        xm = (x1 + x2) / 2
        ax_a.annotate("", xy=(x2 - 0.008, arrow_y), xytext=(x1 + 0.008, arrow_y), arrowprops=dict(arrowstyle="->", lw=2.1, color=palette["line"]))
        ax_a.text(
            xm,
            arrow_y + label_tiers[i],
            arrow_labels[i],
            ha="center",
            va="bottom",
            fontsize=fs_arrow if not compact else fs_arrow * 0.96,
            color=palette["line"],
            bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="#CBD5E1", alpha=0.98),
        )
    direction_text = "Direction lock: Steiger→forward" if compact else "Direction lock: Steiger supports forward"
    ax_a.text(0.05, 0.30 if not compact else 0.29, direction_text, fontsize=fs_small, color="#374151")
    ko_y = 0.185 if not compact else 0.125
    ko_h = 0.16 if not compact else 0.145
    ko_box = patches.FancyBboxPatch((0.245, ko_y), 0.24, ko_h, boxstyle="round,pad=0.01,rounding_size=0.015", linewidth=1.9, edgecolor=palette["accent"], facecolor="#FFF1F2")
    ax_a.add_patch(ko_box)
    ax_a.text(0.365, ko_y + ko_h * 0.69, "Perturbable node", ha="center", va="center", fontsize=fs_arrow, fontweight="bold", color=palette["accent"])
    ax_a.text(0.365, ko_y + ko_h * 0.30, "KO weakens hotspot structure", ha="center", va="center", fontsize=fs_small, color=palette["accent"])
    ax_a.annotate("", xy=(0.53, node_y), xytext=(0.41, 0.34 if not compact else 0.25), arrowprops=dict(arrowstyle="->", lw=1.9, linestyle="--", color=palette["accent"]))
    ax_b.text(0.01, 0.93, "Panel B  Spatial evidence", fontsize=fs_panel, fontweight="bold")
    rect_c1q = patches.FancyBboxPatch((0.08, 0.52), 0.32, 0.3, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=1.8, edgecolor="#7C2D12", facecolor=palette["c1q"], alpha=0.90)
    rect_spp1 = patches.FancyBboxPatch((0.56, 0.52), 0.32, 0.3, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=1.8, edgecolor="#1E3A8A", facecolor=palette["spp1"], alpha=0.90)
    rect_cd8 = patches.FancyBboxPatch((0.33, 0.16), 0.34, 0.20, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=1.8, edgecolor="#065F46", facecolor=palette["cd8"], alpha=0.94)
    ax_b.add_patch(rect_c1q)
    ax_b.add_patch(rect_spp1)
    ax_b.add_patch(rect_cd8)
    ax_b.text(0.24, 0.67, "C1Q hotspot", ha="center", va="center", fontsize=fs_layer, color="white", fontweight="bold")
    ax_b.text(0.72, 0.67, "SPP1 exclusion niche", ha="center", va="center", fontsize=fs_layer, color="white", fontweight="bold")
    ax_b.text(0.50, 0.26, "CD8 interaction zone", ha="center", va="center", fontsize=fs_layer * 0.96, color="white", fontweight="bold")
    ax_b.annotate(
        "attraction ↑",
        xy=(0.40, 0.52),
        xytext=(0.24, 0.40),
        fontsize=fs_arrow,
        color="#7C2D12",
        arrowprops=dict(arrowstyle="->", lw=2, color="#7C2D12"),
    )
    ax_b.annotate(
        "repulsion ↓",
        xy=(0.58, 0.50),
        xytext=(0.69, 0.40),
        fontsize=fs_arrow,
        color="#1E3A8A",
        arrowprops=dict(arrowstyle="-|>", lw=2, color="#1E3A8A"),
    )
    bottom_line = "KO effect on Moran's I varies by dataset (approx. Δ = -0.10 to -0.25)\nscale-dependent interaction: 50/100/200 μm" if compact else "KO effect on Moran's I varies by dataset (approx. Δ = -0.10 to -0.25)  |  scale-dependent interaction: 50/100/200 μm"
    ax_b.text(0.02, 0.03, bottom_line, fontsize=fs_small, color="#111827")
    ax_c.text(0.01, 0.93, "Panel C  Evidence layers", fontsize=fs_panel, fontweight="bold")
    layers = [
        ("Genetic", "MR + bidirectional reverse MR + Steiger"),
        ("Spatial", "Visium + CosMx, Moran's I, niche antagonism"),
        ("Functional", "Virtual KO, dose-response, mediation"),
        ("Clinical", "Survival stratification + ICB association"),
    ]
    layer_h = 0.145 if compact else 0.155
    layer_step = 0.177 if compact else 0.188
    y = 0.75 if not compact else 0.745
    for name, desc in layers:
        r = patches.FancyBboxPatch((0.05, y), 0.90, layer_h, boxstyle="round,pad=0.012,rounding_size=0.015", linewidth=1.4, edgecolor="#334155", facecolor=palette["neutral"])
        ax_c.add_patch(r)
        ax_c.text(0.10, y + layer_h * 0.72, name, fontsize=fs_layer, fontweight="bold", color="#0F172A")
        ax_c.text(0.10, y + layer_h * 0.25, desc, fontsize=fs_arrow, color="#1E293B")
        y -= layer_step
    ax_c.annotate("", xy=(0.5, 0.12), xytext=(0.5, 0.20), arrowprops=dict(arrowstyle="->", lw=1.9, color="#0F172A"))
    ax_c.text(
        0.5,
        0.06,
        "Computationally perturbable axis:\nspatial organization → clinical outcome",
        ha="center",
        va="center",
        fontsize=fs_arrow,
        fontweight="bold",
        color="#0F172A",
    )
    fig.subplots_adjust(top=0.92, bottom=0.06, left=0.03, right=0.98)
    saved_png = safe_savefig(fig, out_png, dpi=350, bbox_inches="tight")
    saved_pdf = safe_savefig(fig, out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)
    return saved_png, saved_pdf


def main():
    root = Path(r"e:\ZhouFX")
    out_dir = ensure_dir(root / "投稿文件" / "main_figures_code_rebuild_from_original")
    _, out1 = draw_evidence_driven_figure1(
        out_dir / "Figure1_rebuilt.png",
        out_dir / "Figure1_rebuilt.pdf",
        figsize=(16, 10),
        scale=1.0,
    )
    _, out2 = draw_evidence_driven_figure1(
        out_dir / "Figure1_rebuilt_2col.png",
        out_dir / "Figure1_rebuilt_2col.pdf",
        figsize=(7.2, 5.1),
        scale=0.63,
    )
    print(out1)
    print(out2)


if __name__ == "__main__":
    main()
