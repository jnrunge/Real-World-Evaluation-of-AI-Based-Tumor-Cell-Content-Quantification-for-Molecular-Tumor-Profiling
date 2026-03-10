from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from math import ceil
from typing import List

import fitz  # PyMuPDF
from PIL import Image
from docx import Document
from docx.shared import Cm


# -----------------------------
# User configuration (edit this)
# -----------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PDF_BASE_DIR = SCRIPT_DIR / "output"           # PDFs relative to this script
OUTPUT_DOCX = SCRIPT_DIR / "output/SI_output.docx"  # Output relative to this script
TEMP_IMAGE_DIR = SCRIPT_DIR / "_si_tmp"      # Temporary image files


@dataclass
class PdfPart:
    pdf: str               # e.g. "set1/plot_A.pdf"
    pages: str = "all"     # e.g. "all", "1", "1-3,5"
    label: str = ""        # optional panel label


@dataclass
class FigureSpec:
    title: str
    description: str
    parts: List[PdfPart]
    cols: int = 1          # layout columns for combined panel image
    width_cm: float = 16.0 # width in Word


FIGURES: List[FigureSpec] = [
    FigureSpec(
        title="Figure S1: Univariate plots for Patho vs FMI TCC discrepancy",
        description="Overview of relationship between [\"Patho\" - \"FMI\"] TCC discrepancy and all investigated variables.",
        parts=[
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_FMI_part1.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_FMI_part2.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_FMI_part3.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_FMI_part4.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_FMI_part5.pdf")
        ],
        cols=1,
        width_cm=12.0,
    ),
    FigureSpec(
        title="Figure S2: Univariate plots for FMI vs AI TCC discrepancy",
        description="Overview of relationship between [\"FMI\" - \"AI\"] TCC discrepancy and all investigated variables.",
        parts=[
            PdfPart(pdf="univar/all_vs_TCC_FMI_minus_TCC_AI_part1.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_FMI_minus_TCC_AI_part2.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_FMI_minus_TCC_AI_part3.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_FMI_minus_TCC_AI_part4.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_FMI_minus_TCC_AI_part5.pdf")
        ],
        cols=1,
        width_cm=12.0,
    ),
    FigureSpec(
        title="Figure S3: Univariate plots for Patho vs AI TCC discrepancy",
        description="Overview of relationship between [\"Patho\" - \"AI\"] TCC discrepancy and all investigated variables.",
        parts=[
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_AI_part1.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_AI_part2.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_AI_part3.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_AI_part4.pdf"),
            PdfPart(pdf="univar/all_vs_TCC_Patho_minus_TCC_AI_part5.pdf")
        ],
        cols=1,
        width_cm=12.0,
    ),
    FigureSpec(
        title="Figure S4: Stability of variable inclusion in bootstrapped models",
        description="Stability of variable inclusion in final models throughout 1,000 bootstrapped repetitions of the model selection process.",
        parts=[
            PdfPart(pdf="bootstrapped_models/bootstrapped_models_variable_inclusion_TCC_FMI_minus_TCC_AI.pdf"),
            PdfPart(pdf="bootstrapped_models/bootstrapped_models_variable_inclusion_TCC_Patho_minus_TCC_AI.pdf"),
            PdfPart(pdf="bootstrapped_models/bootstrapped_models_variable_inclusion_TCC_Patho_minus_TCC_FMI.pdf")
        ],
        cols=1,
        width_cm=12.0,
    ),
    FigureSpec(
        title="Figure S5: Predictor presence in final models",
        description="Overview of the 32 \"final models\" from the model selection process. Multiple models are possible due to randomness in model selection process. More than one final model indicates less certainty about predictor stability.",
        parts=[
            PdfPart(pdf="model_plots/predictor_presence_TCC_Patho_minus_TCC_FMI.pdf"),
            PdfPart(pdf="model_plots/predictor_presence_TCC_FMI_minus_TCC_AI.pdf"),
            PdfPart(pdf="model_plots/predictor_presence_TCC_Patho_minus_TCC_AI.pdf")
        ],
        cols=1,
        width_cm=12.0,
    )
    # Add more FigureSpec(...) entries here
]


# -----------------------------
# Core logic
# -----------------------------
def parse_pages(page_expr: str, total_pages: int) -> List[int]:
    expr = page_expr.strip().lower()
    if expr == "all":
        return list(range(total_pages))

    pages: List[int] = []
    for token in expr.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a, b = token.split("-", 1)
            start = int(a)
            end = int(b)
            if start > end:
                start, end = end, start
            pages.extend(range(start - 1, end))
        else:
            pages.append(int(token) - 1)

    valid = sorted(set(p for p in pages if 0 <= p < total_pages))
    return valid


def render_pdf_pages(pdf_path: Path, page_expr: str, dpi: int = 220) -> List[Image.Image]:
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF not found: {pdf_path}")

    doc = fitz.open(pdf_path)
    try:
        page_indices = parse_pages(page_expr, doc.page_count)
        scale = dpi / 72.0
        matrix = fitz.Matrix(scale, scale)

        images: List[Image.Image] = []
        for i in page_indices:
            pix = doc[i].get_pixmap(matrix=matrix, alpha=False)
            img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
            images.append(img)
        return images
    finally:
        doc.close()


def combine_images(images: List[Image.Image], cols: int = 1, padding: int = 24) -> Image.Image:
    if not images:
        raise ValueError("No images to combine.")

    cols = max(1, cols)
    rows = ceil(len(images) / cols)

    cell_w = max(img.width for img in images)
    cell_h = max(img.height for img in images)

    canvas_w = cols * cell_w + (cols + 1) * padding
    canvas_h = rows * cell_h + (rows + 1) * padding
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

    for idx, img in enumerate(images):
        r = idx // cols
        c = idx % cols
        x0 = padding + c * (cell_w + padding)
        y0 = padding + r * (cell_h + padding)
        x = x0 + (cell_w - img.width) // 2
        y = y0 + (cell_h - img.height) // 2
        canvas.paste(img, (x, y))

    return canvas


def build_docx(figures: List[FigureSpec]) -> None:
    TEMP_IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DOCX.parent.mkdir(parents=True, exist_ok=True)

    doc = Document()
    doc.add_heading("Supplementary Information", level=1)

    for idx, fig in enumerate(figures, start=1):
        doc.add_paragraph(fig.title).runs[0].bold = True

        inserted_any = False
        image_counter = 1

        for part in fig.parts:
            pdf_path = PDF_BASE_DIR / part.pdf
            images = render_pdf_pages(pdf_path, part.pages)

            for img in images:
                inserted_any = True
                panel_file = TEMP_IMAGE_DIR / f"figure_{idx:03d}_{image_counter:03d}.png"
                img.save(panel_file, format="PNG")
                doc.add_picture(str(panel_file), width=Cm(fig.width_cm))
                image_counter += 1

        if inserted_any:
            doc.add_paragraph(fig.description)

    doc.save(OUTPUT_DOCX)


def main() -> None:
    build_docx(FIGURES)
    print(f"Done: {OUTPUT_DOCX}")


if __name__ == "__main__":
    main()
