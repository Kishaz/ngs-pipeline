#!/usr/bin/env python3
"""
Ancestry QC Report
Combines superpopulation proportions table + barplot figure into a PDF.
Designed to be called via Snakemake `script:` directive.
"""

import csv
import os
import glob
from datetime import datetime


def build_combined_tsv(ancestry_dir, tsv_path):
    """Merge per-sample superpop TSVs into one combined file."""
    superpop_files = sorted(
        glob.glob(os.path.join(ancestry_dir, "*.ancestry_superpops.tsv"))
    )
    if not superpop_files:
        raise FileNotFoundError(
            f"No *.ancestry_superpops.tsv files found in {ancestry_dir}"
        )

    # Read header from first file
    with open(superpop_files[0]) as f:
        header = f.readline().strip()

    with open(tsv_path, "w") as out:
        out.write(header + "\n")
        for sp_file in superpop_files:
            with open(sp_file) as f:
                next(f)  # skip header
                for line in f:
                    if line.strip():
                        out.write(line)

    return superpop_files


def generate_pdf(tsv_path, pdf_path, project_name, logo_path=None, barplot_path=None):
    """Generate a 2-page landscape PDF: superpop table + barplot."""
    with open(tsv_path) as f:
        rows = list(csv.reader(f, delimiter="\t"))
    header = rows[0]
    data = rows[1:]

    try:
        from fpdf import FPDF
    except ImportError:
        txt_path = pdf_path.replace(".pdf", ".txt")
        pops = header[1:]
        with open(txt_path, "w") as out:
            out.write("=" * 70 + "\n")
            out.write(f"{project_name} - Ancestry Report\n")
            out.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
            out.write("=" * 70 + "\n\n")
            for row in data:
                vals = [float(v) for v in row[1:]]
                dominant = pops[vals.index(max(vals))]
                out.write(
                    f"{row[0]:20s}  Dominant: {dominant} "
                    f"({max(vals) * 100:.1f}%)\n"
                )
                for p, v in zip(pops, vals):
                    out.write(f"  {p}: {float(v) * 100:.2f}%\n")
                out.write("\n")
        print(f"fpdf2 not installed. Text report: {txt_path}")
        with open(pdf_path, "wb") as f:
            f.write(b"")
        return

    pdf = FPDF(orientation="L", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)

    # ===================== PAGE 1: TABLE =====================
    pdf.add_page()

    if logo_path and os.path.isfile(logo_path):
        pdf.image(logo_path, x=10, y=8, w=20)

    pdf.set_font("Helvetica", "B", 18)
    pdf.cell(
        0, 12, f"{project_name} - Ancestry Report",
        new_x="LMARGIN", new_y="NEXT", align="C",
    )
    pdf.set_font("Helvetica", "", 10)
    pdf.cell(
        0, 6,
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}    |    "
        f"Reference panel: 1000 Genomes    |    Samples: {len(data)}",
        new_x="LMARGIN", new_y="NEXT", align="C",
    )
    pdf.ln(8)

    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(
        0, 8, "Superpopulation Admixture Proportions",
        new_x="LMARGIN", new_y="NEXT",
    )
    pdf.ln(2)

    # Table layout
    pops = header[1:]
    col_w_sample = 40
    col_w_pop = 30
    col_w_dom = 45
    col_ws = [col_w_sample] + [col_w_pop] * len(pops) + [col_w_dom]
    row_h = 8

    # Header
    labels = ["Sample"] + pops + ["Dominant"]
    pdf.set_font("Helvetica", "B", 10)
    pdf.set_fill_color(41, 65, 122)
    pdf.set_text_color(255, 255, 255)
    for i, label in enumerate(labels):
        pdf.cell(col_ws[i], row_h, label, border=1, fill=True, align="C")
    pdf.ln()

    # Rows
    pdf.set_font("Helvetica", "", 9)
    for idx, row in enumerate(data):
        if idx % 2 == 0:
            pdf.set_fill_color(240, 240, 240)
        else:
            pdf.set_fill_color(255, 255, 255)

        vals = [float(v) for v in row[1:]]
        max_idx = vals.index(max(vals))
        dominant = f"{pops[max_idx]} ({vals[max_idx] * 100:.1f}%)"

        # Sample name
        pdf.set_text_color(0, 0, 0)
        pdf.cell(col_ws[0], row_h, row[0], border=1, fill=True, align="C")

        # Population proportions
        for i, v in enumerate(vals):
            pct = f"{v * 100:.2f}%"
            if i == max_idx:
                pdf.set_font("Helvetica", "B", 9)
                pdf.set_text_color(41, 65, 122)
            else:
                pdf.set_font("Helvetica", "", 9)
                pdf.set_text_color(0, 0, 0)
            pdf.cell(col_ws[1 + i], row_h, pct, border=1, fill=True, align="C")

        # Dominant column
        pdf.set_font("Helvetica", "B", 9)
        pdf.set_text_color(41, 65, 122)
        pdf.cell(col_ws[-1], row_h, dominant, border=1, fill=True, align="C")
        pdf.ln()

    # Footer
    pdf.set_text_color(0, 0, 0)
    pdf.ln(4)
    pdf.set_font("Helvetica", "I", 8)
    pdf.cell(
        0, 5,
        "Proportions estimated via ADMIXTURE (supervised, K=5). "
        "Dominant ancestry highlighted in bold.",
        new_x="LMARGIN", new_y="NEXT",
    )
    pdf.cell(
        0, 5,
        "Superpopulations: AFR=African, AMR=Admixed American, "
        "EAS=East Asian, EUR=European, SAS=South Asian",
        new_x="LMARGIN", new_y="NEXT",
    )

    # ===================== PAGE 2: BARPLOT =====================
    if barplot_path and os.path.isfile(barplot_path):
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 14)
        pdf.cell(
            0, 10, "Ancestry Composition",
            new_x="LMARGIN", new_y="NEXT", align="C",
        )
        pdf.ln(4)

        # Scale image to fit page while preserving aspect ratio
        try:
            from PIL import Image
            img = Image.open(barplot_path)
            img_w, img_h = img.size
            aspect = img_h / img_w
        except ImportError:
            aspect = 0.6  # reasonable default

        avail_w = pdf.w - 30
        avail_h = pdf.h - pdf.get_y() - 20

        w = avail_w
        h = w * aspect
        if h > avail_h:
            h = avail_h
            w = h / aspect

        x = (pdf.w - w) / 2
        pdf.image(barplot_path, x=x, y=pdf.get_y(), w=w, h=h)

    pdf.output(pdf_path)
    print(f"PDF report written to: {pdf_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    ancestry_dir = snakemake.params.ancestry_dir  # noqa: F821
    pdf_path = snakemake.output.pdf  # noqa: F821
    project_name = snakemake.params.project_name  # noqa: F821
    logo_path = snakemake.params.get("logo", "")  # noqa: F821
    # Use the container-produced summary (declared as Snakemake input)
    tsv_path = snakemake.input.summary  # noqa: F821
    barplot_path = snakemake.input.barplot  # noqa: F821
else:
    import sys
    ancestry_dir = sys.argv[1]
    pdf_path = sys.argv[2]
    project_name = sys.argv[3] if len(sys.argv) > 3 else "NGS Pipeline"
    logo_path = sys.argv[4] if len(sys.argv) > 4 else ""
    tsv_path = os.path.join(ancestry_dir, "ancestry_summary_superpops.tsv")
    barplot_path = os.path.join(ancestry_dir, "ancestry_barplot_combined.png")

# If the container summary is missing or empty, rebuild from per-sample files
if not os.path.isfile(tsv_path) or os.path.getsize(tsv_path) == 0:
    print("Container summary not found — rebuilding from per-sample files")
    build_combined_tsv(ancestry_dir, tsv_path)

print(f"Using ancestry summary: {tsv_path}")

# Generate report
generate_pdf(tsv_path, pdf_path, project_name, logo_path, barplot_path)
