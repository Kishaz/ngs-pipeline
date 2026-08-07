# =============================================================================
# DNA Alignment: bwa-mem2 (WES / WGS samples only)
# bwa-mem2 is a conda env; samtools is a module
# =============================================================================

rule bwa_mem2_align:
    input:
        reads=get_trimmed_fastqs,
        ref=config["bwa_mem2_index"],
        idx=multiext(
            config["bwa_mem2_index"],
            ".amb", ".ann", ".pac", ".bwt.2bit.64", ".0123",
        ),
    output:
        bam=temp("{results_dir}/alignment/bam/{sample}.raw.bam"),
    params:
        extra=config["alignment"]["bwa_mem2"]["extra"],
        rg=lambda wc: (
            f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}"
            f"\\tLB:{get_seq_type(wc.sample).upper()}"
            f"\\tPL:ILLUMINA\\tPU:{wc.sample}"
        ),
        activate=get_activate_cmd("bwa_mem2"),
    wildcard_constraints:
        sample="|".join([s for s in DNA_SAMPLES if s not in BAM_SAMPLES]) or "NONE",
    envmodules:
        *get_tool_modules("samtools"),
    threads: 30
    resources:
        mem_mb=48000,
        runtime=480,
    log:
        "{results_dir}/logs/bwa_mem2/{sample}.log",
    shell:
        """
        set -o pipefail
        {params.activate}

        bwa-mem2 mem \
            -t {threads} \
            -R '{params.rg}' \
            {params.extra} \
            {input.ref} {input.reads} \
            2> {log} \
        | samtools view -@ 2 -bS -o {output.bam} -

        echo "bwa-mem2 alignment complete" >> {log}
        """
