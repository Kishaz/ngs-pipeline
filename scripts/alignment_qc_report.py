#!/usr/bin/env python3
"""
Alignment QC Report
Parses samtools flagstat + GATK MarkDuplicates metrics → TSV + PDF.
Designed to be called via Snakemake `script:` directive.
"""

import csv
import os
import re
import glob
from datetime import datetime


def parse_flagstat(flagstat_path):
    """Extract alignment metrics from a samtools flagstat file."""
    with open(flagstat_path) as f:
        lines = f.readlines()

    total = 0
    mapped = 0
    paired_ok = 0

    for line in lines:
        if "in total" in line:
            total = int(line.split()[0])
        elif re.search(r"mapped \(", line) and "primary" not in line:
            mapped = int(line.split()[0])
        elif "properly paired" in line:
            paired_ok = int(line.split()[0])

    return total, mapped, paired_ok


def parse_dup_metrics(dup_path):
    """Extract duplication rate from GATK MarkDuplicates metrics."""
    if not os.path.isfile(dup_path):
        return "NA"

    in_data = False
    header = []
    with open(dup_path) as f:
        for line in f:
            if line.startswith("## METRICS CLASS"):
                in_data = True
                continue
            if in_data and not header:
                header = line.strip().split("\t")
                continue
            if in_data and header and line.strip() and not line.startswith("#"):
                vals = line.strip().split("\t")
                if len(vals) > 8:
                    try:
                        return f"{float(vals[8]) * 100:.2f}"
                    except (ValueError, IndexError):
                        return "NA"
    return "NA"


def write_tsv(samples, tsv_path):
    """Write alignment QC metrics to TSV."""
    header = [
        "Sample", "Total_Reads", "Mapped_Reads",
        "Mapping_Rate(%)", "Properly_Paired(%)", "Duplicate_Rate(%)",
    ]
    with open(tsv_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for s in samples:
            w.writerow([
                s["sample"], s["total"], s["mapped"],
                s["map_rate"], s["pair_rate"], s["dup_rate"],
            ])


def generate_pdf(samples, pdf_path, project_name, logo_path=None):
    """Generate a landscape PDF report with alignment metrics table."""
    header = [
        "Sample", "Total Reads", "Mapped Reads",
        "Mapping Rate(%)", "Properly Paired(%)", "Duplicate Rate(%)",
    ]

    try:
        from fpdf import FPDF
    except ImportError:
        txt_path = pdf_path.replace(".pdf", ".txt")
        with open(txt_path, "w") as out:
            out.write("=" * 70 + "\n")
            out.write(f"{project_name} - Alignment QC Report\n")
            out.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
            out.write(f"Reference: GRCh38 v48\n")
            out.write("=" * 70 + "\n\n")
            for s in samples:
                out.write(f"--- {s['sample']} ---\n")
                out.write(f"  Total Reads:       {int(s['total']):,}\n")
                out.write(f"  Mapped Reads:      {int(s['mapped']):,}\n")
                flag = "  ** LOW **" if s["map_rate"] != "NA" and float(s["map_rate"]) < 90 else ""
                out.write(f"  Mapping Rate:      {s['map_rate']}%{flag}\n")
                out.write(f"  Properly Paired:   {s['pair_rate']}%\n")
                flag = "  ** HIGH **" if s["dup_rate"] != "NA" and float(s["dup_rate"]) > 30 else ""
                out.write(f"  Duplicate Rate:    {s['dup_rate']}%{flag}\n\n")
        print(f"fpdf2 not installed. Text report: {txt_path}")
        with open(pdf_path, "wb") as f:
            f.write(b"")
        return

    pdf = FPDF(orientation="L", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    if logo_path and os.path.isfile(logo_path):
        pdf.image(logo_path, x=10, y=8, w=20)

    pdf.set_font("Helvetica", "B", 18)
    pdf.cell(
        0, 12, f"{project_name} - Alignment QC Report",
        new_x="LMARGIN", new_y="NEXT", align="C",
    )
    pdf.set_font("Helvetica", "", 10)
    pdf.cell(
        0, 6,
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}    |    "
        f"Reference: GRCh38 v48    |    Samples: {len(samples)}",
        new_x="LMARGIN", new_y="NEXT", align="C",
    )
    pdf.ln(8)

    col_w = [55, 40, 40, 40, 45, 40]
    row_h = 8

    # Header
    pdf.set_font("Helvetica", "B", 10)
    pdf.set_fill_color(41, 65, 122)
    pdf.set_text_color(255, 255, 255)
    for i, h in enumerate(header):
        pdf.cell(col_w[i], row_h, h, border=1, fill=True, align="C")
    pdf.ln()

    # Data
    pdf.set_font("Helvetica", "", 9)
    for idx, s in enumerate(samples):
        if idx % 2 == 0:
            pdf.set_fill_color(240, 240, 240)
        else:
            pdf.set_fill_color(255, 255, 255)

        row = [
            s["sample"],
            f"{int(s['total']):,}",
            f"{int(s['mapped']):,}",
            s["map_rate"],
            s["pair_rate"],
            s["dup_rate"],
        ]

        for i, val in enumerate(row):
            color_red = False
            if i == 3 and val != "NA" and float(val) < 90:
                color_red = True
            if i == 5 and val != "NA" and float(val) > 30:
                color_red = True

            pdf.set_text_color(200, 30, 30) if color_red else pdf.set_text_color(0, 0, 0)
            pdf.cell(col_w[i], row_h, val, border=1, fill=True, align="C")
        pdf.ln()

    pdf.set_text_color(0, 0, 0)
    pdf.ln(6)
    pdf.set_font("Helvetica", "I", 8)
    pdf.cell(
        0, 5,
        "Flags: Mapping rate < 90% or duplicate rate > 30% highlighted in red.",
        new_x="LMARGIN", new_y="NEXT",
    )
    pdf.cell(
        0, 5,
        "Aligners: bwa-mem2 (DNA) / HISAT2 (RNA)  |  Duplicates: GATK MarkDuplicates",
        new_x="LMARGIN", new_y="NEXT",
    )

    pdf.output(pdf_path)
    print(f"PDF report written to: {pdf_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    metrics_dir = snakemake.params.metrics_dir  # noqa: F821
    tsv_path = snakemake.output.tsv  # noqa: F821
    pdf_path = snakemake.output.pdf  # noqa: F821
    project_name = snakemake.params.project_name  # noqa: F821
    logo_path = snakemake.params.get("logo", "")  # noqa: F821
else:
    import sys
    metrics_dir = sys.argv[1]
    tsv_path = sys.argv[2]
    pdf_path = sys.argv[3]
    project_name = sys.argv[4] if len(sys.argv) > 4 else "NGS Pipeline"
    logo_path = sys.argv[5] if len(sys.argv) > 5 else ""

# Collect metrics
flagstat_files = sorted(glob.glob(os.path.join(metrics_dir, "*.flagstat.txt")))
if not flagstat_files:
    raise FileNotFoundError(f"No flagstat files found in {metrics_dir}")

samples_data = []
for fs_path in flagstat_files:
    sample = os.path.basename(fs_path).replace(".flagstat.txt", "")
    total, mapped, paired_ok = parse_flagstat(fs_path)

    dup_path = os.path.join(metrics_dir, f"{sample}.dup_metrics.txt")
    dup_rate = parse_dup_metrics(dup_path)

    map_rate = f"{mapped / total * 100:.2f}" if total > 0 else "0.00"
    pair_rate = f"{paired_ok / total * 100:.2f}" if total > 0 else "0.00"

    samples_data.append({
        "sample": sample,
        "total": total,
        "mapped": mapped,
        "map_rate": map_rate,
        "pair_rate": pair_rate,
        "dup_rate": dup_rate,
    })
    print(f"  Parsed: {sample}")

write_tsv(samples_data, tsv_path)
print(f"TSV summary written to: {tsv_path}")

generate_pdf(samples_data, pdf_path, project_name, logo_path)
