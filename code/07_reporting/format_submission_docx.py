from pathlib import Path
import argparse

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.shared import Cm, Pt


def apply_style(doc: Document) -> None:
    for section in doc.sections:
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin = Cm(2.5)
        section.right_margin = Cm(2.5)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    normal.paragraph_format.line_spacing = 1.15
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.space_before = Pt(0)

    for style_name, size in [("Heading 1", 16), ("Heading 2", 14), ("Heading 3", 12)]:
        style = doc.styles[style_name]
        style.font.name = "Times New Roman"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
        style.font.size = Pt(size)
        style.font.bold = True
        style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
        style.paragraph_format.line_spacing = 1.15
        style.paragraph_format.space_before = Pt(12)
        style.paragraph_format.space_after = Pt(6)

    lead_labels = ("Background:", "Methods:", "Results:", "Conclusions:", "Keywords:")

    for i, p in enumerate(doc.paragraphs):
        if i == 0 or p.style.name == "Title":
            for r in p.runs:
                r.font.name = "Times New Roman"
                r._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
                r.font.size = Pt(16)
                r.font.bold = True
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.space_after = Pt(18)
            p.paragraph_format.space_before = Pt(0)
            continue

        if p.text.strip() == "Abstract":
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.space_before = Pt(6)
            p.paragraph_format.space_after = Pt(10)
            continue

        if p.text.strip().startswith(lead_labels):
            for j, r in enumerate(p.runs):
                r.font.name = "Times New Roman"
                r._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
                if j == 0:
                    r.font.bold = True
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            p.paragraph_format.line_spacing = 1.15
            if p.text.strip().startswith("Keywords:"):
                p.paragraph_format.space_before = Pt(8)
                p.paragraph_format.space_after = Pt(12)
            else:
                p.paragraph_format.space_after = Pt(6)
            continue

        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
        p.paragraph_format.line_spacing = 1.15
        p.paragraph_format.space_after = Pt(6)


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply submission-style formatting to a DOCX manuscript.")
    parser.add_argument("doc_path", type=str, help="Path to the DOCX file to format.")
    args = parser.parse_args()

    doc_path = Path(args.doc_path).resolve()
    doc = Document(doc_path)
    apply_style(doc)
    doc.save(doc_path)


if __name__ == "__main__":
    main()
