# =============================================================================
# Unified NGS Pipeline
# Supports: RNAseq, WES, WGS | Paired-End & Single-End | FASTQ & BAM input
#
# DAG (DNA):  QC → Align (bwa-mem2) → Sort → MarkDup → BQSR → HaplotypeCaller
#             → GenomicsDBImport → GenotypeGVCFs → Ancestry → Reports
#
# DAG (RNA):  QC → Align (HISAT2)  → Sort → MarkDup → featureCounts
#             → Ancestry → Reports
#
# BAM input:  Auto-detected stage → enters pipeline at the appropriate step
# =============================================================================

import pandas as pd
import os
import subprocess
from snakemake.utils import min_version

min_version("8.0")

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Auto-discover samples from input_dir (if configured)
# ---------------------------------------------------------------------------
_input_dir = (config.get("input_dir") or "").strip()
if _input_dir:
    _discovery_script = os.path.join(workflow.basedir, "scripts", "generate_samples.py")
    _seq_type = config.get("default_seq_type", "rnaseq")
    _samples_tsv = config["samples"]
    subprocess.run(
        [
            sys.executable, _discovery_script,
            "--input-dir", _input_dir,
            "--seq-type", _seq_type,
            "--output", _samples_tsv,
        ],
        check=True,
    )

# ---------------------------------------------------------------------------
# Load sample manifest
# ---------------------------------------------------------------------------
import sys

samples_df = (
    pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    .fillna("")
    .query("sample != ''")
)

# Backward-compatible: add optional columns if missing
if "bam" not in samples_df.columns:
    samples_df["bam"] = ""
if "bam_stage" not in samples_df.columns:
    samples_df["bam_stage"] = ""

SAMPLES = samples_df["sample"].tolist()

# Partition by input type
BAM_SAMPLES = samples_df[
    samples_df["bam"].str.strip() != ""
]["sample"].tolist()

FASTQ_SAMPLES = [s for s in SAMPLES if s not in BAM_SAMPLES]

# Partition by sequencing type
DNA_SAMPLES = samples_df[
    samples_df["seq_type"].isin(["wes", "wgs"])
]["sample"].tolist()

RNA_SAMPLES = samples_df[
    samples_df["seq_type"] == "rnaseq"
]["sample"].tolist()

# PE/SE only applies to FASTQ samples (BAM samples have no R1/R2)
PE_SAMPLES = samples_df[
    (samples_df["R2"].str.strip() != "") & (samples_df["bam"].str.strip() == "")
]["sample"].tolist()

SE_SAMPLES = [s for s in FASTQ_SAMPLES if s not in PE_SAMPLES]

# ---------------------------------------------------------------------------
# Include rule modules (order matters for ruleorder directives)
# ---------------------------------------------------------------------------
include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/align_dna.smk"
include: "rules/align_rna.smk"
include: "rules/post_align.smk"
include: "rules/bqsr.smk"
include: "rules/variant_calling.smk"
include: "rules/quantification.smk"
include: "rules/bam_input.smk"
include: "rules/ancestry.smk"
include: "rules/reports.smk"

# ---------------------------------------------------------------------------
# Staged target rules
# ---------------------------------------------------------------------------
# Each stage depends on the previous stage completing fully (all samples).
# This prevents aggregation rules from running before all inputs exist.
#
# Usage:
#   snakemake --profile profiles/local                  # full pipeline
#   snakemake --profile profiles/local --until stage_qc # QC only
#   snakemake --profile profiles/local --until stage_align  # through alignment
# ---------------------------------------------------------------------------
rd = config["results_dir"]


# --- Stage 1: Read QC (FASTQ samples only) ---
rule stage_qc:
    input:
        [f"{rd}/qc/multiqc/multiqc_report.html",
         f"{rd}/qc/reports/fastq_qc_report.pdf"] if FASTQ_SAMPLES else [],


# --- Stage 2: Alignment + post-alignment + alignment report ---
rule stage_align:
    input:
        rules.stage_qc.input,
        expand(f"{rd}/alignment/bam/{{s}}.sorted.markdup.bam.bai", s=SAMPLES),
        expand(f"{rd}/alignment/metrics/{{s}}.flagstat.txt", s=SAMPLES),
        f"{rd}/alignment/reports/alignment_qc_report.pdf",


# --- Stage 3: BQSR (DNA only; no-op for RNA-only cohorts) ---
rule stage_bqsr:
    input:
        rules.stage_align.input,
        expand(f"{rd}/alignment/bam/{{s}}.sorted.markdup.recal.bam.bai", s=DNA_SAMPLES)
            if DNA_SAMPLES and bqsr_enabled() else [],


# --- Stage 4: Analysis (variant calling / quantification) ---
rule stage_analysis:
    input:
        rules.stage_bqsr.input,
        [f"{rd}/variants/joint/cohort.vcf.gz"] if DNA_SAMPLES else [],
        [f"{rd}/counts/gene_counts_matrix.tsv"] if RNA_SAMPLES else [],


# --- Stage 5: Ancestry estimation ---
rule stage_ancestry:
    input:
        rules.stage_analysis.input,
        [f"{rd}/ancestry/ancestry_summary_superpops.tsv",
         f"{rd}/ancestry/ancestry_barplot_combined.png",
         f"{rd}/ancestry/reports/ancestry_report.pdf"]
            if config["ancestry"]["enabled"] else [],


# ---------------------------------------------------------------------------
# Master rule — chains all stages
# ---------------------------------------------------------------------------
rule all:
    input:
        rules.stage_ancestry.input,


# ---------------------------------------------------------------------------
# Pipeline info
# ---------------------------------------------------------------------------
onstart:
    print("=" * 60)
    print("  NGS Pipeline")
    print("=" * 60)
    print(f"  Samples:     {len(SAMPLES)}")
    print(f"    DNA (WES/WGS): {len(DNA_SAMPLES)}  {DNA_SAMPLES}")
    print(f"    RNA:           {len(RNA_SAMPLES)}  {RNA_SAMPLES}")
    if FASTQ_SAMPLES:
        print(f"    FASTQ input:   {len(FASTQ_SAMPLES)}  (PE: {len(PE_SAMPLES)}, SE: {len(SE_SAMPLES)})")
    if BAM_SAMPLES:
        print(f"    BAM input:     {len(BAM_SAMPLES)}  {BAM_SAMPLES}")
    print(f"  BQSR:        {'enabled' if bqsr_enabled() else 'disabled'}")
    print(f"  Ancestry:    {'enabled' if config['ancestry']['enabled'] else 'disabled'}")
    print(f"  Output:      {rd}")
    print("=" * 60)


onsuccess:
    print("\n" + "=" * 60)
    print("  ALL STAGES COMPLETED SUCCESSFULLY")
    print("=" * 60)
    if FASTQ_SAMPLES:
        print(f"  QC Report:        {rd}/qc/reports/fastq_qc_report.pdf")
    print(f"  Alignment Report: {rd}/alignment/reports/alignment_qc_report.pdf")
    if DNA_SAMPLES:
        if bqsr_enabled():
            print(f"  BQSR:             applied to {len(DNA_SAMPLES)} DNA sample(s)")
        print(f"  Joint VCF:        {rd}/variants/joint/cohort.vcf.gz")
    if RNA_SAMPLES:
        print(f"  Counts Matrix:    {rd}/counts/gene_counts_matrix.tsv")
    if config["ancestry"]["enabled"]:
        print(f"  Ancestry Report:  {rd}/ancestry/reports/ancestry_report.pdf")
    print("=" * 60)


onerror:
    print("\n" + "=" * 60)
    print(f"  PIPELINE FAILED — check logs in {rd}/logs/")
    print(f"  Re-run with:")
    print(f"    SLURM:  snakemake --profile profiles/slurm --rerun-incomplete")
    print(f"    Local:  snakemake --profile profiles/local")
    print("=" * 60)
