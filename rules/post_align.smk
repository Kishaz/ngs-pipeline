# =============================================================================
# Post-Alignment Processing (shared by DNA and RNA)
# samtools sort → GATK MarkDuplicates → samtools index → flagstat
# =============================================================================


rule samtools_sort:
    input:
        "{results_dir}/alignment/bam/{sample}.raw.bam",
    output:
        temp("{results_dir}/alignment/bam/{sample}.sorted.bam"),
    params:
        mem=config["alignment"]["samtools_sort"]["mem_per_thread"],
    envmodules:
        *get_tool_modules("samtools"),
    threads: 30
    resources:
        mem_mb=48000,
        runtime=180,
    log:
        "{results_dir}/logs/samtools_sort/{sample}.log",
    shell:
        """
        samtools sort \
            -@ {threads} \
            -m {params.mem} \
            -o {output} \
            {input} \
            2> {log}
        """


rule mark_duplicates:
    input:
        "{results_dir}/alignment/bam/{sample}.sorted.bam",
    output:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        metrics="{results_dir}/alignment/metrics/{sample}.dup_metrics.txt",
    params:
        extra=config["alignment"]["mark_duplicates"]["extra"],
    envmodules:
        *get_tool_modules("gatk"),
    threads: 1
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        "{results_dir}/logs/markdup/{sample}.log",
    shell:
        """
        gatk MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX false \
            {params.extra} \
            2> {log}
        """


rule samtools_index:
    input:
        "{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
    output:
        "{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
    envmodules:
        *get_tool_modules("samtools"),
    threads: 8
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        """
        samtools index -@ {threads} {input}
        """


rule samtools_flagstat:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
    output:
        "{results_dir}/alignment/metrics/{sample}.flagstat.txt",
    envmodules:
        *get_tool_modules("samtools"),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    shell:
        """
        samtools flagstat {input.bam} > {output}
        """
