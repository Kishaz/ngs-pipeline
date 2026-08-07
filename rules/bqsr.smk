# =============================================================================
# Base Quality Score Recalibration — BQSR (DNA samples only)
#
# 1. BaseRecalibrator  → recalibration model from known variant sites
# 2. ApplyBQSR         → corrected base qualities → analysis-ready BAM
#
# Only runs when bqsr.enabled=true AND known_sites are configured.
# RNA samples ALWAYS skip BQSR.
# =============================================================================


rule base_recalibrator:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
        ref=config["ref"]["fasta"],
        known_sites=config.get("known_sites", []),
    output:
        table="{results_dir}/alignment/bqsr/{sample}.recal_data.table",
    params:
        known_sites_args=lambda wc, input: " ".join(
            f"--known-sites {ks}" for ks in input.known_sites
        ),
        intervals=(
            f"-L {config['variant_calling']['haplotypecaller']['intervals']}"
            if config["variant_calling"]["haplotypecaller"]["intervals"]
            else ""
        ),
    wildcard_constraints:
        sample="|".join(DNA_SAMPLES) if DNA_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("gatk"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=360,
    log:
        "{results_dir}/logs/bqsr/{sample}.base_recalibrator.log",
    shell:
        """
        gatk BaseRecalibrator \
            -R {input.ref} \
            -I {input.bam} \
            {params.known_sites_args} \
            {params.intervals} \
            -O {output.table} \
            2> {log}
        """


rule apply_bqsr:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
        ref=config["ref"]["fasta"],
        table="{results_dir}/alignment/bqsr/{sample}.recal_data.table",
    output:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.recal.bam",
    wildcard_constraints:
        sample="|".join(DNA_SAMPLES) if DNA_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("gatk"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=360,
    log:
        "{results_dir}/logs/bqsr/{sample}.apply_bqsr.log",
    shell:
        """
        gatk ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} \
            2> {log}
        """


rule index_recal_bam:
    input:
        "{results_dir}/alignment/bam/{sample}.sorted.markdup.recal.bam",
    output:
        "{results_dir}/alignment/bam/{sample}.sorted.markdup.recal.bam.bai",
    wildcard_constraints:
        sample="|".join(DNA_SAMPLES) if DNA_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("samtools"),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        """
        samtools index -@ {threads} {input}
        """


# NOTE: No separate flagstat for recalibrated BAMs.
# BQSR only adjusts base quality scores — it does not change read positions,
# mapping status, or duplicate flags. The markdup flagstat is therefore
# identical and is used for the alignment QC report.
