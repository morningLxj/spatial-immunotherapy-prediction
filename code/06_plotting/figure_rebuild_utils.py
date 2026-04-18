from pathlib import Path

import fitz
import matplotlib.pyplot as plt


def render_pdf_first_page(pdf_path, dpi=300):
    doc = fitz.open(str(pdf_path))
    page = doc.load_page(0)
    scale = dpi / 72
    mat = fitz.Matrix(scale, scale)
    pix = page.get_pixmap(matrix=mat, alpha=False)
    img = pix.tobytes("png")
    doc.close()
    return img


def save_single_panel_figure(image_png_bytes, title, out_png, out_pdf, figsize=(12, 8)):
    from PIL import Image
    import io
    import numpy as np

    arr = np.array(Image.open(io.BytesIO(image_png_bytes)))
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(arr)
    ax.axis("off")
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)
    fig.tight_layout()
    fig.savefig(out_png, dpi=350, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=350, bbox_inches="tight")
    plt.close(fig)


def ensure_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
