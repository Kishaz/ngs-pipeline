# =============================================================================
# Common helper functions for sample lookup, SE/PE detection, BAM routing,
# and module/conda tool activation
# =============================================================================

import pandas as pd
import os


# =============================================================================
# TOOL ACTIVATION HELPERS
# =============================================================================
# The config defines each tool as either "module" or "conda".
# - module: loaded via Snakemake envmodules (module load <name>)
# - conda:  activated in the shell block (source activate <name>)
#
# get_tool_modules()  → returns list of module names for envmodules: directive
# get_activate_cmd()  → returns shell string to activate conda envs
# =============================================================================

def get_tool_modules(*tool_keys):
    """
    Return a list of module names for tools that are type=module.
    Conda tools are excluded (they are activated in the shell block).
    Usage in rule: envmodules: *get_tool_modules("fastp", "samtools")
    """
    mods = []
    for key in tool_keys:
        tool = config["tools"].get(key, {})
        if tool.get("type") == "module" and tool.get("name"):
            mods.append(tool["name"])
    return mods


def get_activate_cmd(*tool_keys):
    """
    Return a shell command string that activates any conda environments
    needed for the given tool keys. Module tools are skipped (handled by envmodules).
    Returns empty string if no conda activation needed.
    """
    cmds = []
    needs_anaconda = False
    for key in tool_keys:
        tool = config["tools"].get(key, {})
        if tool.get("type") == "conda" and tool.get("name"):
            needs_anaconda = True
            cmds.append(f"source activate {tool['name']}")
    if not cmds:
        return ""
    # Ensure anaconda module is loaded first (provides conda/source activate)
    prefix = ""
    anaconda = config["tools"].get("anaconda3", {})
    if anaconda.get("type") == "module" and anaconda.get("name"):
        prefix = f"module load {anaconda['name']} && "
    return prefix + " && ".join(cmds) + " && "


# =============================================================================
# SAMPLE LOOKUP FUNCTIONS
# =============================================================================

def get_sample_row(sample):
    """Return the row from samples_df for a given sample name."""
    return samples_df.loc[samples_df["sample"] == sample].iloc[0]


def is_paired(sample):
    """Check if a sample has R2 reads (paired-end)."""
    row = get_sample_row(sample)
    return str(row["R2"]).strip() != ""


def get_seq_type(sample):
    """Return the seq_type for a sample (rnaseq, wes, wgs)."""
    return get_sample_row(sample)["seq_type"]


def is_dna(sample):
    """True if seq_type is WES or WGS."""
    return get_seq_type(sample) in ("wes", "wgs")


def is_rna(sample):
    """True if seq_type is rnaseq."""
    return get_seq_type(sample) == "rnaseq"


def bqsr_enabled():
    """Check if BQSR is enabled in config."""
    return config.get("bqsr", {}).get("enabled", False) and len(config.get("known_sites", [])) > 0


# =============================================================================
# FASTQ INPUT FUNCTIONS
# =============================================================================

def get_raw_r1(wildcards):
    """Return raw R1 path for a sample."""
    return get_sample_row(wildcards.sample)["R1"]


def get_raw_r2(wildcards):
    """Return raw R2 path for a sample."""
    return get_sample_row(wildcards.sample)["R2"]


def get_trimmed_fastqs(wildcards):
    """Return list of trimmed FASTQ paths for a sample."""
    rd = config["results_dir"]
    r1 = f"{rd}/qc/fastp/{wildcards.sample}_R1.trimmed.fastq.gz"
    if is_paired(wildcards.sample):
        return [r1, f"{rd}/qc/fastp/{wildcards.sample}_R2.trimmed.fastq.gz"]
    return [r1]


def get_trimmed_r1(wildcards):
    """Return trimmed R1 path."""
    return f"{config['results_dir']}/qc/fastp/{wildcards.sample}_R1.trimmed.fastq.gz"


def get_trimmed_r2(wildcards):
    """Return trimmed R2 path (only for PE samples)."""
    return f"{config['results_dir']}/qc/fastp/{wildcards.sample}_R2.trimmed.fastq.gz"


# =============================================================================
# BAM ROUTING FUNCTIONS
# =============================================================================

def get_final_bam(wildcards):
    """
    Return the final analysis-ready BAM for a sample.
    DNA + BQSR enabled → recalibrated BAM
    RNA or BQSR disabled → markdup BAM
    """
    rd = config["results_dir"]
    s = wildcards.sample
    if is_dna(s) and bqsr_enabled():
        return f"{rd}/alignment/bam/{s}.sorted.markdup.recal.bam"
    return f"{rd}/alignment/bam/{s}.sorted.markdup.bam"


def get_final_bai(wildcards):
    """Return the index for the final analysis-ready BAM."""
    return get_final_bam(wildcards) + ".bai"


def get_all_final_bams():
    """Return list of all final BAM paths."""
    rd = config["results_dir"]
    bams = []
    for s in SAMPLES:
        if is_dna(s) and bqsr_enabled():
            bams.append(f"{rd}/alignment/bam/{s}.sorted.markdup.recal.bam")
        else:
            bams.append(f"{rd}/alignment/bam/{s}.sorted.markdup.bam")
    return bams


def get_all_final_bais():
    """Return list of all final BAI paths."""
    return [b + ".bai" for b in get_all_final_bams()]


# =============================================================================
# BAM INPUT HELPERS
# =============================================================================

def is_bam_input(sample):
    """True if this sample enters the pipeline via a pre-existing BAM."""
    return str(get_sample_row(sample).get("bam", "")).strip() != ""


def get_bam_stage(sample):
    """Return the BAM processing stage (raw, sorted, markdup, recal) or None."""
    stage = str(get_sample_row(sample).get("bam_stage", "")).strip()
    return stage if stage else None


def get_input_bam(sample):
    """Return the absolute path of the user-provided BAM."""
    return get_sample_row(sample)["bam"]


def _bam_samples_at_stage(stage):
    """Return list of sample names whose bam_stage matches the given stage."""
    return samples_df[
        (samples_df["bam"].str.strip() != "")
        & (samples_df["bam_stage"].str.strip() == stage)
    ]["sample"].tolist()


# =============================================================================
# ANCESTRY HELPERS
# =============================================================================

def get_ancestry_seq_type():
    """
    Determine the seq_type flag for the ancestry container.
    If all samples share the same seq_type, use the mapped value.
    For mixed cohorts, default to 'exome' (most conservative).
    """
    seq_type_map = config.get("ancestry", {}).get("seq_type_map", {})
    sample_types = set(samples_df["seq_type"].tolist())
    if len(sample_types) == 1:
        raw = sample_types.pop()
        return seq_type_map.get(raw, "exome")
    return "exome"


# =============================================================================
# FASTQC OUTPUT FUNCTIONS
# =============================================================================

def get_all_fastqc_raw():
    """Return all expected raw FastQC output files (FASTQ samples only)."""
    out = []
    rd = config["results_dir"]
    for _, row in samples_df.iterrows():
        s = row["sample"]
        if str(row.get("bam", "")).strip():
            continue  # skip BAM-input samples
        out.append(f"{rd}/qc/fastqc_raw/{s}_R1_fastqc.html")
        if str(row["R2"]).strip():
            out.append(f"{rd}/qc/fastqc_raw/{s}_R2_fastqc.html")
    return out


def get_all_fastqc_trimmed():
    """Return all expected trimmed FastQC output files (FASTQ samples only)."""
    out = []
    rd = config["results_dir"]
    for _, row in samples_df.iterrows():
        s = row["sample"]
        if str(row.get("bam", "")).strip():
            continue  # skip BAM-input samples
        out.append(f"{rd}/qc/fastqc_trimmed/{s}_R1.trimmed_fastqc.html")
        if str(row["R2"]).strip():
            out.append(f"{rd}/qc/fastqc_trimmed/{s}_R2.trimmed_fastqc.html")
    return out
