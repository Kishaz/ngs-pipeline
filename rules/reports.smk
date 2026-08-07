# =============================================================================
# QC Report Generation: FastQ QC, Alignment QC, Ancestry
# =============================================================================


rule fastq_qc_report:
    input:
        jsons=expand(
            "{rd}/qc/fastp/{sample}_fastp.json",
            rd=config["results_dir"],
            sample=FASTQ_SAMPLES,
        ),
    output:
        tsv="{results_dir}/qc/reports/fastq_qc_summary.tsv",
        pdf="{results_dir}/qc/reports/fastq_qc_report.pdf",
    params:
        fastp_dir=lambda wc: f"{wc.results_dir}/qc/fastp",
        project_name=config["reports"]["project_name"],
        logo=config["reports"]["lab_logo"],
    envmodules:
        *get_tool_modules("python3"),
    resources:
        mem_mb=4000,
        runtime=15,
    log:
        "{results_dir}/logs/reports/fastq_qc_report.log",
    script:
        "../scripts/fastq_qc_report.py"


rule alignment_qc_report:
    input:
        flagstats=expand(
            "{rd}/alignment/metrics/{sample}.flagstat.txt",
            rd=config["results_dir"],
            sample=SAMPLES,
        ),
        dup_metrics=expand(
            "{rd}/alignment/metrics/{sample}.dup_metrics.txt",
            rd=config["results_dir"],
            sample=SAMPLES,
        ),
    output:
        tsv="{results_dir}/alignment/reports/alignment_qc_summary.tsv",
        pdf="{results_dir}/alignment/reports/alignment_qc_report.pdf",
    params:
        metrics_dir=lambda wc: f"{wc.results_dir}/alignment/metrics",
        project_name=config["reports"]["project_name"],
        logo=config["reports"]["lab_logo"],
    envmodules:
        *get_tool_modules("python3"),
    resources:
        mem_mb=4000,
        runtime=15,
    log:
        "{results_dir}/logs/reports/alignment_qc_report.log",
    script:
        "../scripts/alignment_qc_report.py"


rule ancestry_qc_report:
    input:
        summary="{results_dir}/ancestry/ancestry_summary_superpops.tsv",
        barplot="{results_dir}/ancestry/ancestry_barplot_combined.png",
    output:
        pdf="{results_dir}/ancestry/reports/ancestry_report.pdf",
    params:
        ancestry_dir=lambda wc: f"{wc.results_dir}/ancestry",
        project_name=config["reports"]["project_name"],
        logo=config["reports"]["lab_logo"],
    envmodules:
        *get_tool_modules("python3"),
    resources:
        mem_mb=4000,
        runtime=15,
    log:
        "{results_dir}/logs/reports/ancestry_report.log",
    script:
        "../scripts/ancestry_qc_report.py"
