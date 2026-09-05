# =============================================================================
# RNA-seq–specific QC (RNA samples only)
#   RSeQC infer_experiment.py  -> library strandedness (drives featureCounts -s)
#   RSeQC read_distribution.py -> read coverage over CDS/UTR/intron/intergenic
#   Picard CollectRnaSeqMetrics -> 5'->3' bias, %coding/%mRNA/%ribosomal, etc.
#
# All three feed MultiQC. The inferred strand is written to a tiny ".strand"
# file (one integer: featureCounts -s value) that featureCounts and Picard read.
#
# RSeQC is a conda env (config tools.rseqc); Picard is a module (tools.picard).
# BED12 / refFlat / rRNA interval_list are prebuilt from the GENCODE GTF by
# scripts/gtf_to_rnaseq_refs.py (see setup_references.sh) — no UCSC tools needed.
# =============================================================================

RNAQC_SCRIPTS = os.path.join(workflow.basedir, "scripts")

_RNA_ONLY = "|".join(RNA_SAMPLES) if RNA_SAMPLES else "NONE"


# ---------------------------------------------------------------------------
# RSeQC — infer strandedness, and write the featureCounts -s value
# ---------------------------------------------------------------------------
rule rseqc_infer_experiment:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
        bed=config.get("rnaseq_qc", {}).get("bed12", ""),
    output:
        txt="{results_dir}/qc/rseqc/{sample}.infer_experiment.txt",
        strand="{results_dir}/qc/rseqc/{sample}.strand",
    params:
        activate=get_activate_cmd("rseqc"),
        parser=f"{RNAQC_SCRIPTS}/infer_strandedness.py",
        sample_size=config.get("rnaseq_qc", {}).get("infer_sample_size", 400000),
        threshold=config.get("rnaseq_qc", {}).get("strand_threshold", 0.80),
    wildcard_constraints:
        sample=_RNA_ONLY,
    threads: 1
    resources:
        mem_mb=4000,
        runtime=45,
    log:
        "{results_dir}/logs/rseqc/{sample}.infer_experiment.log",
    shell:
        """
        {params.activate}
        infer_experiment.py \
            -i {input.bam} \
            -r {input.bed} \
            -s {params.sample_size} \
            > {output.txt} 2> {log}

        python3 {params.parser} --threshold {params.threshold} \
            {output.txt} > {output.strand} 2>> {log}
        """


# ---------------------------------------------------------------------------
# RSeQC — read distribution over gene features
# ---------------------------------------------------------------------------
rule rseqc_read_distribution:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
        bed=config.get("rnaseq_qc", {}).get("bed12", ""),
    output:
        "{results_dir}/qc/rseqc/{sample}.read_distribution.txt",
    params:
        activate=get_activate_cmd("rseqc"),
    wildcard_constraints:
        sample=_RNA_ONLY,
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "{results_dir}/logs/rseqc/{sample}.read_distribution.log",
    shell:
        """
        {params.activate}
        read_distribution.py \
            -i {input.bam} \
            -r {input.bed} \
            > {output} 2> {log}
        """


# ---------------------------------------------------------------------------
# Picard — CollectRnaSeqMetrics (5'->3' bias, coding/UTR/intronic/ribosomal)
# STRAND_SPECIFICITY is taken from the inferred strand (or the numeric config):
#   featureCounts -s 0 -> NONE
#                 -s 1 -> SECOND_READ_TRANSCRIPTION_STRAND
#                 -s 2 -> FIRST_READ_TRANSCRIPTION_STRAND
# ---------------------------------------------------------------------------
rule picard_rnaseqmetrics:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
        refflat=config.get("rnaseq_qc", {}).get("refflat", ""),
        strand=get_strand_input,
    output:
        "{results_dir}/qc/picard/{sample}.rna_metrics.txt",
    params:
        strandedness=config["quantification"]["featurecounts"]["strandedness"],
        rrna=config.get("rnaseq_qc", {}).get("rrna_intervals", ""),
    wildcard_constraints:
        sample=_RNA_ONLY,
    envmodules:
        # Picard's CollectRnaSeqMetrics is invoked via the GATK module (bundles
        # Picard + its Java); the standalone picard wrapper on this cluster does
        # not dispatch subcommands reliably. samtools pre-filters the stream.
        *get_tool_modules("gatk", "samtools"),
    threads: 2
    resources:
        # CollectRnaSeqMetrics holds per-transcript coverage; on a large
        # BAM it needs real heap. Pin -Xmx (GATK otherwise sizes to the
        # NODE's RAM, not the SLURM cgroup, and OOMs) and leave headroom
        # under mem_mb for JVM non-heap + the samtools pre-filter.
        mem_mb=16000,
        runtime=120,
    log:
        "{results_dir}/logs/picard/{sample}.rna_metrics.log",
    shell:
        """
        # Resolve featureCounts-style strand (0/1/2): numeric config, else inferred
        S="{params.strandedness}"
        if [ "$S" = "auto" ]; then S=$(cat {input.strand}); fi
        case "$S" in
            1) SPEC=SECOND_READ_TRANSCRIPTION_STRAND ;;
            2) SPEC=FIRST_READ_TRANSCRIPTION_STRAND ;;
            *) SPEC=NONE ;;
        esac

        RRNA=""
        if [ -n "{params.rrna}" ] && [ -f "{params.rrna}" ]; then
            RRNA="--RIBOSOMAL_INTERVALS {params.rrna}"
        fi

        # CollectRnaSeqMetrics computes pair orientation for the stranded
        # metrics and throws on inter-chromosomal read pairs (mate on a
        # different contig) — common in tumour/FFPE RNA-seq. Drop only those
        # reads on the way in: keep unpaired reads, mate-unmapped reads, and
        # pairs whose mate is on the same contig. Stream to gatk via stdin.
        samtools view -h -b \
            -e '!flag.paired || flag.munmap || rnext == "=" || rnext == rname' \
            {input.bam} 2> {log} \
        | gatk --java-options "-Xmx12g" CollectRnaSeqMetrics \
            -I /dev/stdin \
            -O {output} \
            --REF_FLAT {input.refflat} \
            --STRAND_SPECIFICITY $SPEC \
            $RRNA \
            2>> {log}
        """
