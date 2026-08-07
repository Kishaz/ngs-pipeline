# =============================================================================
# BAM Input Ingestion Rules
# Symlink user-provided BAMs into the pipeline's expected output paths,
# allowing samples to enter the pipeline at the appropriate stage.
#
# Stage detection is based on the bam_stage column in samples.tsv:
#   raw     → enters at samtools_sort
#   sorted  → enters at mark_duplicates
#   markdup → enters at indexing (+ flagstat, quantification, variant calling)
#   recal   → enters at indexing (DNA samples with BQSR already applied)
# =============================================================================

import os


# ---------------------------------------------------------------------------
# Helper: list of samples at each BAM stage
# ---------------------------------------------------------------------------
_RAW_BAM_SAMPLES = _bam_samples_at_stage("raw")
_SORTED_BAM_SAMPLES = _bam_samples_at_stage("sorted")
_MARKDUP_BAM_SAMPLES = _bam_samples_at_stage("markdup")
_RECAL_BAM_SAMPLES = _bam_samples_at_stage("recal")


# ---------------------------------------------------------------------------
# Ruleorder: ingestion rules take priority over alignment/processing rules
# for BAM-input samples (wildcard_constraints ensure no overlap)
# ---------------------------------------------------------------------------
ruleorder: ingest_raw_bam > bwa_mem2_align
ruleorder: ingest_raw_bam > hisat2_align
ruleorder: ingest_sorted_bam > samtools_sort
ruleorder: ingest_markdup_bam > mark_duplicates
ruleorder: ingest_recal_bam > apply_bqsr


# ---------------------------------------------------------------------------
# Stage: raw → pipeline enters at samtools_sort
# ---------------------------------------------------------------------------
rule ingest_raw_bam:
    output:
        bam=temp("{results_dir}/alignment/bam/{sample}.raw.bam"),
    wildcard_constraints:
        sample="|".join(_RAW_BAM_SAMPLES) if _RAW_BAM_SAMPLES else "NONE",
    run:
        src = os.path.abspath(get_input_bam(wildcards.sample))
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        if os.path.islink(output.bam):
            os.remove(output.bam)
        os.symlink(src, output.bam)


# ---------------------------------------------------------------------------
# Stage: sorted → pipeline enters at mark_duplicates
# ---------------------------------------------------------------------------
rule ingest_sorted_bam:
    output:
        bam="{results_dir}/alignment/bam/{sample}.sorted.bam",
    wildcard_constraints:
        sample="|".join(_SORTED_BAM_SAMPLES) if _SORTED_BAM_SAMPLES else "NONE",
    run:
        src = os.path.abspath(get_input_bam(wildcards.sample))
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        if os.path.islink(output.bam):
            os.remove(output.bam)
        os.symlink(src, output.bam)


# ---------------------------------------------------------------------------
# Stage: markdup → pipeline enters at indexing + flagstat
# Also creates a placeholder dup_metrics so alignment QC report works.
# ---------------------------------------------------------------------------
rule ingest_markdup_bam:
    output:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        metrics="{results_dir}/alignment/metrics/{sample}.dup_metrics.txt",
    wildcard_constraints:
        sample="|".join(_MARKDUP_BAM_SAMPLES) if _MARKDUP_BAM_SAMPLES else "NONE",
    run:
        src = os.path.abspath(get_input_bam(wildcards.sample))
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        os.makedirs(os.path.dirname(output.metrics), exist_ok=True)
        if os.path.islink(output.bam):
            os.remove(output.bam)
        os.symlink(src, output.bam)

        # Symlink index if available alongside the input BAM
        bai_src = src + ".bai"
        if not os.path.isfile(bai_src):
            bai_src = src.replace(".bam", ".bai")
        if os.path.isfile(bai_src):
            bai_dst = output.bam + ".bai"
            if os.path.islink(bai_dst):
                os.remove(bai_dst)
            os.symlink(os.path.abspath(bai_src), bai_dst)

        # Placeholder dup_metrics (Picard header-only format)
        with open(output.metrics, "w") as f:
            f.write("## BAM provided pre-processed (markdup stage)\n")
            f.write("## No duplicate marking metrics available\n")
            f.write("LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
                    "SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t"
                    "UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
                    "READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t"
                    "ESTIMATED_LIBRARY_SIZE\n")
            f.write(f"{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")


# ---------------------------------------------------------------------------
# Stage: recal → pipeline enters at indexing (DNA + BQSR already applied)
# Also creates placeholder markdup BAM symlink + dup_metrics.
# ---------------------------------------------------------------------------
rule ingest_recal_bam:
    output:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.recal.bam",
        markdup_bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        metrics="{results_dir}/alignment/metrics/{sample}.dup_metrics.txt",
    wildcard_constraints:
        sample="|".join(_RECAL_BAM_SAMPLES) if _RECAL_BAM_SAMPLES else "NONE",
    run:
        src = os.path.abspath(get_input_bam(wildcards.sample))
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        os.makedirs(os.path.dirname(output.metrics), exist_ok=True)

        # Symlink recal BAM
        if os.path.islink(output.bam):
            os.remove(output.bam)
        os.symlink(src, output.bam)

        # Symlink index if available
        bai_src = src + ".bai"
        if not os.path.isfile(bai_src):
            bai_src = src.replace(".bam", ".bai")
        if os.path.isfile(bai_src):
            bai_dst = output.bam + ".bai"
            if os.path.islink(bai_dst):
                os.remove(bai_dst)
            os.symlink(os.path.abspath(bai_src), bai_dst)

        # Also create markdup BAM symlink (same file — for flagstat/reports)
        if os.path.islink(output.markdup_bam):
            os.remove(output.markdup_bam)
        os.symlink(src, output.markdup_bam)

        # Placeholder dup_metrics
        with open(output.metrics, "w") as f:
            f.write("## BAM provided pre-processed (recal stage)\n")
            f.write("## No duplicate marking metrics available\n")
            f.write("LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
                    "SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t"
                    "UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
                    "READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t"
                    "ESTIMATED_LIBRARY_SIZE\n")
            f.write(f"{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")
