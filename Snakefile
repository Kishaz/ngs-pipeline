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

# ---------------------------------------------------------------------------
# Resolve results_dir.
#
#   - results_dir ABSOLUTE  -> used as-is.
#   - results_dir RELATIVE  -> placed UNDER results_base (so results_base is the
#                              parent for named runs, e.g. base=/scratch/me +
#                              results_dir=run1 -> /scratch/me/run1).
#   - results_dir EMPTY     -> auto-generate results_<timestamp> under
#                              results_base, PERSISTED to .results_dir so resumes
#                              (which re-parse this file) reuse the SAME folder
#                              instead of minting a new timestamp (which would
#                              silently restart the pipeline). Delete .results_dir
#                              to start a genuinely new auto-dated run.
#
#   results_base defaults to the ngs_pipeline repo directory; set it to a
#   scratch path so relative/auto-dated results land there.
# ---------------------------------------------------------------------------
_results_dir = (config.get("results_dir") or "").strip()
_base = (config.get("results_base") or workflow.basedir).rstrip("/")
if _results_dir:
    if not os.path.isabs(_results_dir):
        _results_dir = os.path.join(_base, _results_dir)
else:
    _marker = os.path.join(workflow.basedir, ".results_dir")
    if os.path.isfile(_marker):
        with open(_marker) as _f:
            _results_dir = _f.read().strip()
    if not _results_dir:
        from datetime import datetime
        _results_dir = os.path.join(_base, f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        with open(_marker, "w") as _f:
            _f.write(_results_dir)
    print(f"[results_dir] not provided — using {_results_dir} (persisted in .results_dir)")
config["results_dir"] = _results_dir

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
include: "rules/rnaseq_qc.smk"
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
# MultiQC + RNA-specific QC (RSeQC/Picard) live here so they run for every
# cohort — including BAM-input RNA, where stage_qc (FASTQ read QC) is a no-op.
rule stage_align:
    input:
        rules.stage_qc.input,
        expand(f"{rd}/alignment/bam/{{s}}.sorted.markdup.bam.bai", s=SAMPLES),
        expand(f"{rd}/alignment/metrics/{{s}}.flagstat.txt", s=SAMPLES),
        f"{rd}/alignment/reports/alignment_qc_report.pdf",
        [f"{rd}/qc/multiqc/multiqc_report.html"] if SAMPLES else [],
        get_rnaseq_qc_targets(),


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
    default_target: True
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

    # --- Reproducibility: record pipeline git commit + resolved references ---
    try:
        import datetime as _dt
        def _git(*a):
            try:
                return subprocess.check_output(
                    ["git", "-C", workflow.basedir, *a],
                    stderr=subprocess.DEVNULL,
                ).decode().strip()
            except Exception:
                return "unknown"
        os.makedirs(rd, exist_ok=True)
        _rq = config.get("rnaseq_qc", {})
        with open(os.path.join(rd, "provenance.txt"), "w") as _p:
            _p.write("NGS Pipeline — run provenance\n")
            _p.write(f"timestamp:        {_dt.datetime.now().isoformat(timespec='seconds')}\n")
            _p.write(f"pipeline_commit:  {_git('rev-parse', 'HEAD')}\n")
            _p.write(f"pipeline_branch:  {_git('rev-parse', '--abbrev-ref', 'HEAD')}\n")
            _p.write(f"pipeline_dirty:   {'yes' if _git('status', '--porcelain') else 'no'}\n")
            _p.write(f"pipeline_dir:     {workflow.basedir}\n")
            _p.write(f"samples:          {len(SAMPLES)} (DNA {len(DNA_SAMPLES)}, RNA {len(RNA_SAMPLES)})\n")
            _p.write(f"results_dir:      {rd}\n")
            _p.write(f"ref_fasta:        {config['ref']['fasta']}\n")
            _p.write(f"hisat2_index:     {config.get('hisat2_index', '')}\n")
            _p.write(f"gtf:              {config['gtf']}\n")
            _p.write(f"strandedness:     {config['quantification']['featurecounts']['strandedness']}\n")
            _p.write(f"rnaseq_qc:        {'enabled' if _rq.get('enabled', True) else 'disabled'}\n")
            _p.write(f"rnaseq_qc.bed12:  {_rq.get('bed12', '')}\n")
            _p.write(f"rnaseq_qc.refflat:{_rq.get('refflat', '')}\n")
        print(f"  Provenance:  {rd}/provenance.txt")
    except Exception as _e:
        print(f"  WARNING: could not write provenance.txt: {_e}")
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

    # --- Publish small final deliverables off scratch to a persistent dir ---
    # Big intermediates (BAMs, temp files) stay in results_dir (scratch); the
    # small shareable outputs are COPIED to publish_dir so they survive a scratch
    # purge. Defaults to <repo>/deliverables/<run-name>. Set publish_dir: "" to
    # disable, or point it at any persistent path.
    import os as _os, shutil as _sh, glob as _gl
    _pub_cfg = config.get("publish_dir", "__DEFAULT__")
    _pub = (_os.path.join(workflow.basedir, "deliverables",
                          _os.path.basename(rd.rstrip("/")))
            if _pub_cfg == "__DEFAULT__" else (_pub_cfg or "").strip())
    if _pub:
        _items = [
            "provenance.txt",
            "counts/gene_counts_matrix.tsv",
            "qc/reports", "qc/multiqc/multiqc_report.html",
            "qc/rseqc", "qc/picard", "alignment/reports",
            "ancestry/ancestry_summary_superpops.tsv", "ancestry/reports",
            "variants/joint",  # DNA joint VCF (may be large)
        ]
        try:
            for _rel in _items:
                _src, _dst = _os.path.join(rd, _rel), _os.path.join(_pub, _rel)
                if _os.path.isdir(_src):
                    _sh.copytree(_src, _dst, dirs_exist_ok=True)
                elif _os.path.isfile(_src):
                    _os.makedirs(_os.path.dirname(_dst), exist_ok=True)
                    _sh.copy2(_src, _dst)
            for _pat in ("ancestry/*.png", "ancestry/*.svg"):
                for _src in _gl.glob(_os.path.join(rd, _pat)):
                    _dst = _os.path.join(_pub, "ancestry", _os.path.basename(_src))
                    _os.makedirs(_os.path.dirname(_dst), exist_ok=True)
                    _sh.copy2(_src, _dst)
            print(f"  Deliverables:     published to {_pub}")
        except Exception as _e:
            print(f"  WARNING: could not publish deliverables to {_pub}: {_e}")
    print("=" * 60)


onerror:
    print("\n" + "=" * 60)
    print(f"  PIPELINE FAILED — check logs in {rd}/logs/")
    print(f"  Re-run with:")
    print(f"    SLURM:  snakemake --profile profiles/slurm --rerun-incomplete")
    print(f"    Local:  snakemake --profile profiles/local")
    print("=" * 60)
