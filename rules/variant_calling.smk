# =============================================================================
# Variant Calling (DNA samples only)
# GATK HaplotypeCaller → GenomicsDBImport → GenotypeGVCFs
#
# Uses recalibrated BAMs (post-BQSR) when enabled, else markdup BAMs.
# =============================================================================


rule haplotype_caller:
    input:
        bam=get_final_bam,
        bai=get_final_bai,
        ref=config["ref"]["fasta"],
    output:
        gvcf="{results_dir}/variants/gvcf/{sample}.g.vcf.gz",
        tbi="{results_dir}/variants/gvcf/{sample}.g.vcf.gz.tbi",
    params:
        extra=config["variant_calling"]["haplotypecaller"]["extra"],
        intervals=(
            f"-L {config['variant_calling']['haplotypecaller']['intervals']}"
            if config["variant_calling"]["haplotypecaller"]["intervals"]
            else ""
        ),
    wildcard_constraints:
        sample="|".join(DNA_SAMPLES) if DNA_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("gatk"),
    threads: 4
    resources:
        mem_mb=16000,
        runtime=720,
    log:
        "{results_dir}/logs/haplotypecaller/{sample}.log",
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            {params.intervals} \
            {params.extra} \
            --native-pair-hmm-threads {threads} \
            2> {log}
        """


rule genomics_db_import:
    input:
        gvcfs=expand(
            "{rd}/variants/gvcf/{sample}.g.vcf.gz",
            rd=config["results_dir"],
            sample=DNA_SAMPLES,
        ),
    output:
        db=directory("{results_dir}/variants/genomicsdb"),
    params:
        gvcf_args=lambda wc, input: " ".join(f"-V {g}" for g in input.gvcfs),
        intervals=config["variant_calling"]["genomicsdb"]["intervals"],
        batch_size=config["variant_calling"]["genomicsdb"]["batch_size"],
    envmodules:
        *get_tool_modules("gatk"),
    threads: 4
    resources:
        mem_mb=48000,
        runtime=480,
    log:
        "{results_dir}/logs/genomicsdb_import.log",
    shell:
        """
        gatk GenomicsDBImport \
            {params.gvcf_args} \
            --genomicsdb-workspace-path {output.db} \
            -L {params.intervals} \
            --batch-size {params.batch_size} \
            --reader-threads {threads} \
            2> {log}
        """


rule genotype_gvcfs:
    input:
        db="{results_dir}/variants/genomicsdb",
        ref=config["ref"]["fasta"],
    output:
        vcf="{results_dir}/variants/joint/cohort.vcf.gz",
        tbi="{results_dir}/variants/joint/cohort.vcf.gz.tbi",
    params:
        extra=config["variant_calling"]["genotypegvcfs"]["extra"],
    envmodules:
        *get_tool_modules("gatk"),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=480,
    log:
        "{results_dir}/logs/genotypegvcfs.log",
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            {params.extra} \
            2> {log}

        gatk IndexFeatureFile -I {output.vcf} 2>> {log}
        """
