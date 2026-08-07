# =============================================================================
# Read QC: fastp (PE/SE), FastQC (raw + trimmed), MultiQC
# =============================================================================

SCRIPTS = os.path.join(workflow.basedir, "scripts")


# ---------------------------------------------------------------------------
# fastp — Paired-End
# ---------------------------------------------------------------------------
rule fastp_pe:
    input:
        R1=get_raw_r1,
        R2=get_raw_r2,
    output:
        r1=temp("{results_dir}/qc/fastp/{sample}_R1.trimmed.fastq.gz"),
        r2=temp("{results_dir}/qc/fastp/{sample}_R2.trimmed.fastq.gz"),
        json="{results_dir}/qc/fastp/{sample}_fastp.json",
        html=temp("{results_dir}/qc/fastp/{sample}_fastp.html"),
    params:
        qual=config["qc"]["fastp"]["qualified_quality_phred"],
        minlen=config["qc"]["fastp"]["length_required"],
        script=f"{SCRIPTS}/run_fastp.sh",
    wildcard_constraints:
        sample="|".join(PE_SAMPLES) if PE_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("fastp"),
    threads: 16
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        "{results_dir}/logs/fastp/{sample}.log",
    shell:
        """
        bash {params.script} \
            --r1 {input.R1} --r2 {input.R2} \
            --out-r1 {output.r1} --out-r2 {output.r2} \
            --json {output.json} --html {output.html} \
            --qual {params.qual} --minlen {params.minlen} \
            --threads {threads} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# fastp — Single-End
# ---------------------------------------------------------------------------
rule fastp_se:
    input:
        R1=get_raw_r1,
    output:
        r1=temp("{results_dir}/qc/fastp/{sample}_R1.trimmed.fastq.gz"),
        json="{results_dir}/qc/fastp/{sample}_fastp.json",
        html=temp("{results_dir}/qc/fastp/{sample}_fastp.html"),
    params:
        qual=config["qc"]["fastp"]["qualified_quality_phred"],
        minlen=config["qc"]["fastp"]["length_required"],
        script=f"{SCRIPTS}/run_fastp.sh",
    wildcard_constraints:
        sample="|".join(SE_SAMPLES) if SE_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("fastp"),
    threads: 16
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        "{results_dir}/logs/fastp/{sample}.log",
    shell:
        """
        bash {params.script} \
            --r1 {input.R1} \
            --out-r1 {output.r1} \
            --json {output.json} --html {output.html} \
            --qual {params.qual} --minlen {params.minlen} \
            --threads {threads} \
            2> {log}
        """


ruleorder: fastp_pe > fastp_se


# ---------------------------------------------------------------------------
# FastQC — Raw reads
# ---------------------------------------------------------------------------
rule fastqc_raw_r1:
    input:
        get_raw_r1,
    output:
        html="{results_dir}/qc/fastqc_raw/{sample}_R1_fastqc.html",
        zip="{results_dir}/qc/fastqc_raw/{sample}_R1_fastqc.zip",
    params:
        outdir=lambda wc: f"{wc.results_dir}/qc/fastqc_raw",
        script=f"{SCRIPTS}/run_fastqc.sh",
    envmodules:
        *get_tool_modules("fastqc"),
    threads: 8
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "{results_dir}/logs/fastqc/{sample}_R1_raw.log",
    shell:
        """
        bash {params.script} \
            --input {input} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}

        base=$(basename {input} .fastq.gz)
        if [ -f "{params.outdir}/${{base}}_fastqc.html" ] && \
           [ "{params.outdir}/${{base}}_fastqc.html" != "{output.html}" ]; then
            mv "{params.outdir}/${{base}}_fastqc.html" "{output.html}"
            mv "{params.outdir}/${{base}}_fastqc.zip" "{output.zip}"
        fi
        """


rule fastqc_raw_r2:
    input:
        get_raw_r2,
    output:
        html="{results_dir}/qc/fastqc_raw/{sample}_R2_fastqc.html",
        zip="{results_dir}/qc/fastqc_raw/{sample}_R2_fastqc.zip",
    params:
        outdir=lambda wc: f"{wc.results_dir}/qc/fastqc_raw",
        script=f"{SCRIPTS}/run_fastqc.sh",
    wildcard_constraints:
        sample="|".join(PE_SAMPLES) if PE_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("fastqc"),
    threads: 8
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "{results_dir}/logs/fastqc/{sample}_R2_raw.log",
    shell:
        """
        bash {params.script} \
            --input {input} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}

        base=$(basename {input} .fastq.gz)
        if [ -f "{params.outdir}/${{base}}_fastqc.html" ] && \
           [ "{params.outdir}/${{base}}_fastqc.html" != "{output.html}" ]; then
            mv "{params.outdir}/${{base}}_fastqc.html" "{output.html}"
            mv "{params.outdir}/${{base}}_fastqc.zip" "{output.zip}"
        fi
        """


# ---------------------------------------------------------------------------
# FastQC — Trimmed reads
# ---------------------------------------------------------------------------
rule fastqc_trimmed_r1:
    input:
        "{results_dir}/qc/fastp/{sample}_R1.trimmed.fastq.gz",
    output:
        html="{results_dir}/qc/fastqc_trimmed/{sample}_R1.trimmed_fastqc.html",
        zip="{results_dir}/qc/fastqc_trimmed/{sample}_R1.trimmed_fastqc.zip",
    params:
        outdir=lambda wc: f"{wc.results_dir}/qc/fastqc_trimmed",
        script=f"{SCRIPTS}/run_fastqc.sh",
    envmodules:
        *get_tool_modules("fastqc"),
    threads: 8
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "{results_dir}/logs/fastqc/{sample}_R1_trimmed.log",
    shell:
        """
        bash {params.script} \
            --input {input} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}
        """


rule fastqc_trimmed_r2:
    input:
        "{results_dir}/qc/fastp/{sample}_R2.trimmed.fastq.gz",
    output:
        html="{results_dir}/qc/fastqc_trimmed/{sample}_R2.trimmed_fastqc.html",
        zip="{results_dir}/qc/fastqc_trimmed/{sample}_R2.trimmed_fastqc.zip",
    params:
        outdir=lambda wc: f"{wc.results_dir}/qc/fastqc_trimmed",
        script=f"{SCRIPTS}/run_fastqc.sh",
    wildcard_constraints:
        sample="|".join(PE_SAMPLES) if PE_SAMPLES else "NONE",
    envmodules:
        *get_tool_modules("fastqc"),
    threads: 8
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "{results_dir}/logs/fastqc/{sample}_R2_trimmed.log",
    shell:
        """
        bash {params.script} \
            --input {input} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# MultiQC
# ---------------------------------------------------------------------------
rule multiqc:
    input:
        fastqc_raw=get_all_fastqc_raw(),
        fastqc_trimmed=get_all_fastqc_trimmed(),
        fastp=expand(
            "{rd}/qc/fastp/{sample}_fastp.json",
            rd=config["results_dir"],
            sample=FASTQ_SAMPLES,
        ),
    output:
        "{results_dir}/qc/multiqc/multiqc_report.html",
    params:
        indir=lambda wc: f"{wc.results_dir}/qc",
        outdir=lambda wc: f"{wc.results_dir}/qc/multiqc",
        script=f"{SCRIPTS}/run_multiqc.sh",
    envmodules:
        *get_tool_modules("python3", "multiqc"),
    resources:
        mem_mb=4000,
        runtime=30,
    log:
        "{results_dir}/logs/multiqc.log",
    shell:
        """
        bash {params.script} \
            --indir {params.indir} \
            --outdir {params.outdir} \
            2> {log}
        """
