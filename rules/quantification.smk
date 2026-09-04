# =============================================================================
# RNA Quantification: featureCounts (RNA samples only)
# featureCounts (subread) is a conda env
# =============================================================================


rule featurecounts:
    input:
        bam="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam",
        bai="{results_dir}/alignment/bam/{sample}.sorted.markdup.bam.bai",
        gtf=config["gtf"],
        strand=get_strand_input,
    output:
        counts="{results_dir}/counts/{sample}.featureCounts.txt",
        summary="{results_dir}/counts/{sample}.featureCounts.txt.summary",
    params:
        extra=config["quantification"]["featurecounts"]["extra"],
        strandedness=config["quantification"]["featurecounts"]["strandedness"],
        activate=get_activate_cmd("subread"),
    wildcard_constraints:
        sample="|".join(RNA_SAMPLES) if RNA_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("samtools"),
    threads: 16
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        "{results_dir}/logs/featurecounts/{sample}.log",
    shell:
        """
        {params.activate}

        # Strandedness: numeric config value used as-is; "auto" reads the value
        # inferred by RSeQC infer_experiment.py (rules/rnaseq_qc.smk).
        S="{params.strandedness}"
        if [ "$S" = "auto" ]; then S=$(cat {input.strand}); fi

        # featureCounts needs -p --countReadPairs for paired-end BAMs and must
        # NOT get it for single-end. The manifest can't tell us for BAM-input
        # samples (no R2), so detect from the BAM's first read flag (bit 0x1).
        FLAG=$(samtools view {input.bam} 2>/dev/null | head -1 | cut -f2 || true)
        PE=""
        if [ -n "$FLAG" ] && [ $(( FLAG & 1 )) -ne 0 ]; then
            PE="-p --countReadPairs"
        fi

        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            {params.extra} \
            -s $S \
            $PE \
            -T {threads} \
            {input.bam} \
            2> {log}
        """


rule merge_counts:
    input:
        expand(
            "{rd}/counts/{sample}.featureCounts.txt",
            rd=config["results_dir"],
            sample=RNA_SAMPLES,
        ),
    output:
        "{results_dir}/counts/gene_counts_matrix.tsv",
    envmodules:
        *get_tool_modules("python3"),
    resources:
        mem_mb=4000,
        runtime=15,
    log:
        "{results_dir}/logs/featurecounts/merge_counts.log",
    run:
        import pandas as pd
        import os

        dfs = []
        for f in input:
            sample_name = os.path.basename(f).replace(".featureCounts.txt", "")
            df = pd.read_csv(f, sep="\t", comment="#")
            df = df[["Geneid", df.columns[-1]]].rename(
                columns={df.columns[-1]: sample_name}
            )
            dfs.append(df)

        merged = dfs[0]
        for df in dfs[1:]:
            merged = merged.merge(df, on="Geneid", how="outer")

        merged = merged.fillna(0).astype(
            {c: int for c in merged.columns if c != "Geneid"}
        )
        merged.to_csv(output[0], sep="\t", index=False)
