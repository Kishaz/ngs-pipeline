# =============================================================================
# Ancestry Estimation
# Stage BAMs → run containerized ancestry pipeline
# Singularity is at /usr/bin/singularity (system binary, no module needed)
# =============================================================================


rule stage_bams_for_ancestry:
    input:
        bam=get_all_final_bams(),
        bai=get_all_final_bais(),
    output:
        staging=directory("{results_dir}/ancestry/.staging"),
        done=touch("{results_dir}/ancestry/.staging_done"),
    run:
        import os, re

        os.makedirs(output.staging, exist_ok=True)
        for bam_path in input.bam:
            # Strip all pipeline suffixes to recover the clean sample name
            basename = os.path.basename(bam_path)
            sample = re.sub(
                r"\.(sorted\.markdup\.recal|sorted\.markdup|sorted|markdup)\.bam$",
                "",
                basename,
            )
            # Sanitize underscores → hyphens (ADMIXTURE requirement)
            safe_name = sample.replace("_", "-")
            # Stage with clean {sample}.bam name so the container sees
            # exactly one entry per sample
            dst_bam = os.path.join(output.staging, f"{safe_name}.bam")
            dst_bai = os.path.join(output.staging, f"{safe_name}.bam.bai")
            src_bam = os.path.abspath(bam_path)
            src_bai = os.path.abspath(bam_path + ".bai")
            if not os.path.exists(dst_bam):
                os.symlink(src_bam, dst_bam)
            if not os.path.exists(dst_bai):
                os.symlink(src_bai, dst_bai)


rule run_ancestry:
    input:
        staging="{results_dir}/ancestry/.staging",
        done="{results_dir}/ancestry/.staging_done",
    output:
        summary="{results_dir}/ancestry/ancestry_summary_superpops.tsv",
    params:
        container=config["ancestry"]["container"],
        seq_type=get_ancestry_seq_type(),
        threads=config["ancestry"]["threads"],
        subpops="--subpops" if config["ancestry"]["subpops"] else "",
        outdir=lambda wc: os.path.abspath(f"{wc.results_dir}/ancestry"),
        staging_abs=lambda wc: os.path.abspath(
            f"{wc.results_dir}/ancestry/.staging"
        ),
        script=os.path.join(workflow.basedir, "scripts", "run_ancestry.sh"),
    # Singularity is at /usr/bin/singularity — no module needed
    threads: config["ancestry"]["threads"]
    resources:
        mem_mb=80000,
        runtime=1440,
    log:
        "{results_dir}/logs/ancestry/ancestry_batch.log",
    shell:
        """
        bash {params.script} \
            --bam-dir {params.staging_abs} \
            --outdir {params.outdir} \
            --container {params.container} \
            --seq-type {params.seq_type} \
            --threads {params.threads} \
            {params.subpops} \
            2>&1 | tee {log}
        """


rule plot_ancestry:
    input:
        summary="{results_dir}/ancestry/ancestry_summary_superpops.tsv",
    output:
        combined_png="{results_dir}/ancestry/ancestry_barplot_combined.png",
        combined_svg="{results_dir}/ancestry/ancestry_barplot_combined.svg",
        superpops_png="{results_dir}/ancestry/ancestry_barplot_superpops.png",
    params:
        ancestry_dir=lambda wc: os.path.abspath(f"{wc.results_dir}/ancestry"),
        script=os.path.join(workflow.basedir, "scripts", "plot_ancestry.R"),
    envmodules:
        *get_tool_modules("R"),
    resources:
        mem_mb=4000,
        runtime=15,
    log:
        "{results_dir}/logs/ancestry/plot_ancestry.log",
    shell:
        """
        Rscript {params.script} {params.ancestry_dir} 2>&1 | tee {log}
        """
