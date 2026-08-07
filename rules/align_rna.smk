# =============================================================================
# RNA Alignment: HISAT2 (RNAseq samples only)
# Both HISAT2 and samtools are modules
# =============================================================================

rule hisat2_align:
    input:
        reads=get_trimmed_fastqs,
        idx=multiext(
            config["hisat2_index"],
            ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
            ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2",
        ),
    output:
        bam=temp("{results_dir}/alignment/bam/{sample}.raw.bam"),
        summary="{results_dir}/alignment/metrics/{sample}.hisat2_summary.txt",
    params:
        idx=config["hisat2_index"],
        extra=config["alignment"]["hisat2"]["extra"],
        input_flags=lambda wc: (
            f"-1 {get_trimmed_r1(wc)} -2 {get_trimmed_r2(wc)}"
            if is_paired(wc.sample)
            else f"-U {get_trimmed_r1(wc)}"
        ),
        rg=lambda wc: (
            f"--rg-id {wc.sample} "
            f"--rg SM:{wc.sample} "
            f"--rg LB:RNASEQ "
            f"--rg PL:ILLUMINA "
            f"--rg PU:{wc.sample}"
        ),
    wildcard_constraints:
        sample="|".join([s for s in RNA_SAMPLES if s not in BAM_SAMPLES]) or "NONE",
    envmodules:
        *get_tool_modules("hisat2", "samtools"),
    threads: 30
    resources:
        mem_mb=32000,
        runtime=480,
    log:
        "{results_dir}/logs/hisat2/{sample}.log",
    shell:
        """
        set -o pipefail

        hisat2 \
            -x {params.idx} \
            {params.input_flags} \
            {params.extra} \
            {params.rg} \
            -p {threads} \
            --new-summary \
            --summary-file {output.summary} \
            2> {log} \
        | samtools view -@ 2 -bS -o {output.bam} -

        echo "HISAT2 alignment complete" >> {log}
        """
