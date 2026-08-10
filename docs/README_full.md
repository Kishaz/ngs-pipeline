# NGS Pipeline — Full Technical Documentation

A production-grade Snakemake pipeline for processing **RNAseq**, **WES**, and **WGS** data through quality control, alignment, base quality recalibration, variant calling / gene quantification, and genetic ancestry estimation.

Supports both **paired-end** and **single-end** reads, **FASTQ** and **pre-aligned BAM** inputs with automatic sample discovery. Each stage writes outputs to dedicated directories in standard formats for downstream analysis.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Prerequisites](#prerequisites)
3. [Reference Resources Setup](#reference-resources-setup)
4. [Configuration](#configuration)
5. [Input Modes](#input-modes)
6. [Pipeline Stages — Detailed](#pipeline-stages--detailed)
7. [Execution Profiles](#execution-profiles)
8. [Re-running Specific Stages](#re-running-specific-stages)
9. [Output Directory Structure](#output-directory-structure)
10. [Storage Optimization](#storage-optimization)
11. [Standalone Script Usage](#standalone-script-usage)
12. [Resource Requirements](#resource-requirements)
13. [Troubleshooting](#troubleshooting)

---

## Pipeline Overview

```
Raw FASTQ / Pre-aligned BAM
  │
  ▼
┌──────────────────────────────────────────────────────┐
│  STAGE 1: READ QC                                    │
│  fastp → FastQC (raw + trimmed) → MultiQC            │
│  Output: trimmed FASTQs (temp), QC reports, PDF      │
└──────────────────────┬───────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────┐
│  STAGE 2: ALIGNMENT                                  │
│  DNA → bwa-mem2       RNA → HISAT2                   │
│  → samtools sort → GATK MarkDuplicates → index       │
│  Output: sorted, deduplicated BAM + flagstat          │
└──────────────────────┬───────────────────────────────┘
                       │
          ┌────────────┴────────────┐
          ▼                         ▼
┌───────────────────┐    ┌────────────────────┐
│  STAGE 3a: BQSR   │    │  (RNA skips BQSR)  │
│  (DNA only)       │    │                    │
│  BaseRecalibrator │    └────────┬───────────┘
│  → ApplyBQSR      │             │
│  Output: recal BAM│             │
└────────┬──────────┘             │
         │                        │
         ▼                        ▼
┌───────────────────┐    ┌────────────────────┐
│  STAGE 4a:        │    │  STAGE 4b:         │
│  VARIANT CALLING  │    │  QUANTIFICATION    │
│  (DNA only)       │    │  (RNA only)        │
│  HaplotypeCaller  │    │  featureCounts     │
│  → GenomicsDB     │    │  → merged matrix   │
│  → GenotypeGVCFs  │    │                    │
│  Output: joint VCF│    │  Output: counts TSV│
└────────┬──────────┘    └────────┬───────────┘
         │                        │
         └────────────┬───────────┘
                      ▼
┌──────────────────────────────────────────────────────┐
│  STAGE 5: ANCESTRY ESTIMATION                        │
│  Stage BAMs → Singularity container (batch mode)     │
│  GATK pileup → ADMIXTURE (K=5/K=23) → PCA           │
│  Output: ancestry proportions, PCA plots, barplots   │
└──────────────────────┬───────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────┐
│  STAGE 6: REPORTS                                    │
│  PDF reports at each stage (fpdf2)                   │
│  Read QC → Alignment QC → Ancestry report            │
└──────────────────────────────────────────────────────┘
```

---

## Prerequisites

### Software Dependencies

| Tool | Version | Stage | Installation |
|------|---------|-------|--------------|
| Snakemake | >= 8.0 | Orchestration | `pip install snakemake` |
| snakemake-executor-plugin-slurm | latest | SLURM job submission | `pip install snakemake-executor-plugin-slurm` |
| fastp | 1.0.1 | Read QC | HPC module or conda |
| FastQC | 0.12.1 | Read QC | HPC module or conda |
| MultiQC | 1.20 | Read QC | HPC module or conda |
| bwa-mem2 | 2.2+ | DNA alignment | Conda environment |
| HISAT2 | 2.2.1 | RNA alignment | HPC module or conda |
| SAMtools | 1.21 | BAM processing | HPC module or conda |
| GATK | 4.6.2.0 | MarkDup, BQSR, variant calling | HPC module or conda |
| Subread (featureCounts) | 2.0+ | RNA quantification | Conda environment |
| Singularity | 3.x | Ancestry container | System installation |
| Python | 3.11+ | Reports, sample discovery | HPC module |
| R | 4.4.2 | Ancestry plots | HPC module (with Bioconductor) |

### Python Packages

```bash
pip install --user fpdf2 pandas Pillow
```

### R Packages

```r
install.packages(c("ggplot2", "gridExtra", "grid"))
```

### Tool Activation Configuration

The pipeline supports two tool activation methods per tool, configured in `config/config.yaml`:

```yaml
tools:
  fastp:
    type: "module"          # "module" for HPC environment modules
    name: "fastp/1.0.1"     # module name to load
  bwa_mem2:
    type: "conda"           # "conda" for conda environments
    name: "bwa-mem2"        # conda environment name
```

**Module-based tools:** fastp, FastQC, MultiQC, HISAT2, SAMtools, GATK, Python, R
**Conda-based tools:** bwa-mem2, Subread (featureCounts)

---

## Reference Resources Setup

All reference files live under a single `resources/` directory. The setup script `scripts/setup_references.sh` automates downloading and indexing.

### Quick Setup

```bash
# Allocate resources for indexing (bwa-mem2 needs ~80 GB RAM, HISAT2 needs ~200 GB)
srun --mem=200G --cpus-per-task=16 --time=06:00:00 --pty bash

# Run setup
bash scripts/setup_references.sh --outdir /path/to/ngs_resources --threads 16
```

### What Gets Downloaded

#### 1. Reference Genome — GRCh38 Primary Assembly

| Item | Detail |
|------|--------|
| **Source** | GENCODE / Ensembl |
| **URL** | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz` |
| **Version** | GRCh38.p14 (GENCODE release 48) |
| **Contents** | Chromosomes 1–22, X, Y, MT, unlocalized/unplaced scaffolds. Excludes ALT contigs and patches. |
| **Size** | ~3.1 GB (uncompressed) |
| **Output** | `resources/genome/GRCh38.primary_assembly.genome.fa` |

#### 2. Gene Annotation — GENCODE v48

| Item | Detail |
|------|--------|
| **Source** | GENCODE |
| **URL** | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz` |
| **Version** | GENCODE v48 (Ensembl 114) |
| **Contents** | Comprehensive gene annotation: protein-coding genes, lncRNAs, pseudogenes, miRNAs. 63,070 genes, 255,824 transcripts. |
| **Size** | ~1.5 GB (uncompressed GTF) |
| **Output** | `resources/annotation/gencode.v48.primary_assembly.annotation.gtf` |

#### 3. Known Variant Sites — GATK Resource Bundle (for BQSR)

These must be downloaded manually from the GATK resource bundle (requires `gsutil`):

| File | Source | Purpose |
|------|--------|---------|
| `Homo_sapiens_assembly38.dbsnp138.vcf.gz` | `gs://genomics-public-data/resources/broad/hg38/v0/` | dbSNP v138 — common human SNPs |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | Same bucket | Validated insertion/deletion sites |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | Same bucket | Additional known indel sites |

**Download command:**
```bash
# Using gsutil (Google Cloud SDK)
gsutil -m cp \
  gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
  gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
  gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
  gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
  /path/to/resources/known_sites/
```

**Alternative web download:**
Browse to `https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0` and download manually.

#### 4. Ancestry Container

The ancestry pipeline runs inside a Singularity container containing GATK 3.8, PLINK 1.9, and ADMIXTURE 1.3.0 with the 1000 Genomes Phase 3 reference panel (2,504 individuals, 26 populations, 5 continental super-populations).

Place the `.sif` file at: `resources/ancestry/ancestry-pipeline.sif`

#### 5. WES Capture Kit Intervals (WES only)

Place your capture kit BED file in `resources/intervals/`. This restricts variant calling and BQSR to on-target regions.

### What Gets Indexed

The setup script builds four types of indices:

| Index | Tool | RAM Required | Time | Output Files |
|-------|------|-------------|------|--------------|
| FASTA index | `samtools faidx` | < 1 GB | < 1 min | `.fa.fai` |
| Sequence dictionary | `gatk CreateSequenceDictionary` | < 4 GB | < 1 min | `.dict` |
| bwa-mem2 index | `bwa-mem2 index` | ~80 GB | 1–2 hours | `.amb`, `.ann`, `.pac`, `.bwt.2bit.64`, `.0123` |
| HISAT2 index | `hisat2-build` | ~200 GB | 2–3 hours | 8 × `.ht2` files + `splice_sites.tsv`, `exons.tsv` |

**HISAT2 splice-aware index:** If a GTF annotation is present, the script extracts splice sites and exons using `hisat2_extract_splice_sites.py` and `hisat2_extract_exons.py`, then builds a splice-aware index. This dramatically improves RNA-seq alignment accuracy at splice junctions.

### Final Resource Directory Structure

```
resources/GRCh38_v48/
├── genome/
│   ├── GRCh38.primary_assembly.genome.fa          Reference FASTA
│   ├── GRCh38.primary_assembly.genome.fa.fai      samtools index
│   ├── GRCh38.primary_assembly.genome.dict        GATK sequence dictionary
│   ├── GRCh38.primary_assembly.genome.fa.amb      ┐
│   ├── GRCh38.primary_assembly.genome.fa.ann      │ bwa-mem2
│   ├── GRCh38.primary_assembly.genome.fa.pac      │ index files
│   ├── GRCh38.primary_assembly.genome.fa.0123     │
│   └── GRCh38.primary_assembly.genome.fa.bwt.2bit.64 ┘
├── hisat2_index/
│   ├── grch38.1.ht2 … grch38.8.ht2               HISAT2 genome index (8 files)
│   ├── splice_sites.tsv                           Extracted from GENCODE GTF
│   └── exons.tsv                                  Extracted from GENCODE GTF
├── annotation/
│   └── gencode.v48.primary_assembly.annotation.gtf
├── known_sites/
│   ├── Homo_sapiens_assembly38.dbsnp138.vcf.gz (+.tbi)
│   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz (+.tbi)
│   └── Homo_sapiens_assembly38.known_indels.vcf.gz (+.tbi)
├── intervals/
│   └── <capture_kit>.bed                          WES target regions (user-provided)
└── ancestry/
    └── ancestry-pipeline.sif                      Singularity container (~2 GB)
```

**Total disk:** ~60 GB

---

## Configuration

### Main Configuration File: `config/config.yaml`

```yaml
# ── Input ──────────────────────────────────────────
input_dir: "/path/to/your/fastqs_or_bams"   # Auto-discovery (empty = manual TSV)
default_seq_type: "rnaseq"                   # rnaseq | wes | wgs
samples: "config/samples.tsv"                # Sample manifest

# ── Reference Paths ────────────────────────────────
ref:
  fasta: "/path/to/resources/genome/GRCh38.primary_assembly.genome.fa"
hisat2_index: "/path/to/resources/hisat2_index/grch38_v48"
bwa_mem2_index: "/path/to/resources/genome/GRCh38.primary_assembly.genome.fa"
gtf: "/path/to/resources/annotation/gencode.v48.primary_assembly.annotation.gtf"

# ── QC Parameters ──────────────────────────────────
qc:
  fastp:
    qualified_quality_phred: 20    # Min base quality to keep
    length_required: 50            # Min read length after trimming

# ── Alignment Parameters ──────────────────────────
alignment:
  bwa_mem2:
    extra: "-Y"                    # Soft-clip supplementary alignments
  hisat2:
    extra: "--dta"                 # Downstream transcript assembly mode
  samtools_sort:
    mem_per_thread: "4G"           # RAM per sorting thread
  mark_duplicates:
    extra: "--REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT"

# ── BQSR (DNA only) ───────────────────────────────
bqsr:
  enabled: true
known_sites:
  - "/path/to/resources/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
  - "/path/to/resources/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  - "/path/to/resources/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"

# ── Variant Calling (DNA only) ────────────────────
variant_calling:
  haplotypecaller:
    extra: ""
    intervals: ""                  # WES: path to capture kit .bed
  genomicsdb:
    intervals: ""                  # Required for joint calling (.list)
    batch_size: 50
  genotypegvcfs:
    extra: ""

# ── RNA Quantification ────────────────────────────
quantification:
  featurecounts:
    extra: "-t exon -g gene_id --primary"
    strandedness: 0                # 0=unstranded, 1=stranded, 2=reverse-stranded

# ── Ancestry ───────────────────────────────────────
ancestry:
  enabled: true
  container: "/path/to/resources/ancestry/ancestry-pipeline.sif"
  runtime: "singularity"           # singularity | docker
  seq_type_map:
    rnaseq: "rna"
    wes: "exome"
    wgs: "wgs"
  subpops: true                    # Enable K=23 sub-population analysis
  threads: 30

# ── Reports ────────────────────────────────────────
reports:
  lab_logo: "resources/Lab_Logo.png"
  project_name: "NGS Pipeline"

# ── Results Directory ──────────────────────────────
results_dir: "results"
```

### Tool Module Configuration

Each tool specifies its activation method (`module` or `conda`) and version:

```yaml
tools:
  fastp:     { type: "module", name: "fastp/1.0.1" }
  fastqc:    { type: "module", name: "fastqc/0.12.1" }
  multiqc:   { type: "module", name: "multiqc/1.20" }
  bwa_mem2:  { type: "conda",  name: "bwa-mem2" }
  hisat2:    { type: "module", name: "HISAT2/2.2.1" }
  samtools:  { type: "module", name: "samtools/1.21" }
  gatk:      { type: "module", name: "gatk/4.6.2.0" }
  subread:   { type: "conda",  name: "subread" }
  python3:   { type: "module", name: "python311" }
  anaconda3: { type: "module", name: "anaconda3/2023.09" }
  R:         { type: "module", name: "R/4.4.2+Bioconductor" }
```

---

## Input Modes

### Option A: Auto-discovery (Recommended)

Set `input_dir` in config and the pipeline auto-discovers samples on each run:

```yaml
input_dir: "/path/to/your/fastqs"
default_seq_type: "rnaseq"
```

**FASTQ detection patterns (priority order):**

| R1 Pattern | R2 Pattern | Example |
|------------|------------|---------|
| `*_R1_001.fastq.gz` | `*_R2_001.fastq.gz` | `SampleA_R1_001.fastq.gz` |
| `*_R1.fastq.gz` | `*_R2.fastq.gz` | `SampleA_R1.fastq.gz` |
| `*_1.fastq.gz` | `*_2.fastq.gz` | `SampleA_1.fastq.gz` |
| `*.R1.fastq.gz` | `*.R2.fastq.gz` | `SampleA.R1.fastq.gz` |

If no matching R2 is found, the sample is classified as single-end.

**BAM stage detection (longest suffix match):**

| Filename Suffix | Detected Stage | Pipeline Entry Point |
|----------------|----------------|---------------------|
| `.sorted.markdup.recal.bam` | `recal` | Index only |
| `.sorted.markdup.bam` | `markdup` | Index + flagstat |
| `.sorted.bam` | `sorted` | MarkDuplicates |
| `.raw.bam` or `.bam` | `raw` | samtools sort |

### Option B: Manual Sample Manifest

Leave `input_dir` empty and create `config/samples.tsv` (tab-delimited):

```
sample	R1	R2	seq_type	bam	bam_stage
PatientA	/data/PatientA_R1.fastq.gz	/data/PatientA_R2.fastq.gz	wes
PatientB	/data/PatientB_R1.fastq.gz		rnaseq
PatientC				wes	/data/PatientC.sorted.markdup.bam	markdup
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `R1` | FASTQ input | Path to R1 FASTQ |
| `R2` | PE only | Path to R2 FASTQ (empty for SE) |
| `seq_type` | Yes | `rnaseq`, `wes`, or `wgs` |
| `bam` | BAM input | Path to pre-aligned BAM (empty for FASTQ) |
| `bam_stage` | BAM input | Processing stage: `raw`, `sorted`, `markdup`, `recal` |

### Standalone Sample Discovery

```bash
python scripts/generate_samples.py \
  --input-dir /path/to/data \
  --seq-type rnaseq \
  --output config/samples.tsv
```

---

## Pipeline Stages — Detailed

### Stage 1: Read Quality Control

#### fastp — Adapter Removal and Quality Trimming

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--qualified_quality_phred` | 20 | Minimum base quality to retain |
| `--length_required` | 50 | Minimum read length after trimming |
| `--detect_adapter_for_pe` | (PE only) | Auto-detect PE adapter sequences |
| Threads | 16 | Parallel processing threads |

**Outputs:** Trimmed FASTQs (temporary — deleted after alignment), JSON report (kept for MultiQC), HTML report (temporary).

**Script:** `scripts/run_fastp.sh` — Detects PE vs SE from presence of `--r2` flag. Validates all inputs exist before execution.

#### FastQC — Quality Assessment

Runs twice per sample:
1. **Pre-trimming:** On raw FASTQ files
2. **Post-trimming:** On fastp-processed files

**Outputs:** HTML reports (kept), ZIP archives (temporary — deleted after MultiQC).

**Script:** `scripts/run_fastqc.sh` — Handles output renaming when FastQC basename differs from expected sample name.

#### MultiQC — Aggregated Report

Combines all FastQC and fastp JSON reports into a single interactive HTML report.

**Script:** `scripts/run_multiqc.sh` — Requires `python311` module loaded before `multiqc/1.20`.

#### Read QC PDF Report

Custom Python script generates a shareable PDF with per-sample QC metrics table.

**Flags (highlighted in red):** Pass rate < 80%, Q30 < 80%, Duplication > 30%.

**Script:** `scripts/fastq_qc_report.py` — Uses fpdf2 with optional lab logo. Falls back to text-only if Pillow not installed.

---

### Stage 2: Sequence Alignment

#### DNA Alignment: bwa-mem2

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-Y` | enabled | Soft-clip supplementary alignments |
| `-t` | 16 | Alignment threads |
| `-R` | `@RG\tID:{sample}\tSM:{sample}\tLB:{SEQ_TYPE}\tPL:ILLUMINA\tPU:{sample}` | Full read group header |

Alignment is piped directly to `samtools view -bS` for immediate BAM conversion.

**Script:** `scripts/run_bwa_mem2.sh` — Validates all bwa-mem2 index files exist (`.amb`, `.ann`, `.pac`, `.bwt.2bit.64` or `.bwt.2bit`). Supports checkpoint files for standalone use.

#### RNA Alignment: HISAT2

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--dta` | enabled | Downstream transcript assembly mode |
| `-p` | 16 | Alignment threads |
| `--new-summary` | enabled | Machine-parseable summary format |
| `--rg-id` | `{sample}` | Read group ID |
| `--rg SM:` | `{sample}` | Sample name in read group |
| `--rg LB:` | `RNASEQ` | Library identifier |
| `--rg PL:` | `ILLUMINA` | Sequencing platform |

**PE mode:** `-1 {R1} -2 {R2}` | **SE mode:** `-U {R1}`

Alignment is piped directly to `samtools view -bS` for immediate BAM conversion.

**Script:** `scripts/run_hisat2.sh` — Detects PE/SE from input flags. Outputs alignment summary alongside BAM.

#### Post-Alignment Processing (all samples)

| Step | Tool | Key Parameters | Output |
|------|------|---------------|--------|
| Coordinate sort | `samtools sort` | `-m 4G` per thread, 30 threads | Sorted BAM (temporary) |
| Mark duplicates | `gatk MarkDuplicates` | `REMOVE_DUPLICATES=false`, `VALIDATION_STRINGENCY=LENIENT` | Markdup BAM + metrics |
| Index | `samtools index` | 8 threads | `.bai` index |
| Stats | `samtools flagstat` | — | Alignment summary text |

---

### Stage 3: Base Quality Score Recalibration (DNA only)

> RNA samples always skip BQSR. BQSR only runs if `bqsr.enabled: true` AND `known_sites` list is non-empty.

| Step | Tool | Input | Output |
|------|------|-------|--------|
| Build model | `gatk BaseRecalibrator` | Markdup BAM + ref + 3 known-sites VCFs | `.recal_data.table` |
| Apply correction | `gatk ApplyBQSR` | Markdup BAM + recal table | Recalibrated BAM |
| Index | `samtools index` | Recalibrated BAM | `.bai` index |

**Known sites files (from GATK resource bundle):**
- `Homo_sapiens_assembly38.dbsnp138.vcf.gz` — Common human SNPs
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` — Validated indels
- `Homo_sapiens_assembly38.known_indels.vcf.gz` — Additional known indels

**WES note:** Set `variant_calling.haplotypecaller.intervals` to your capture kit BED to restrict recalibration to on-target regions.

**Script:** `scripts/run_bqsr.sh` — Validates known sites VCFs exist and are tabix-indexed. Supports optional intervals.

---

### Stage 4a: Variant Calling (DNA only)

| Step | Tool | Key Parameters | Output |
|------|------|---------------|--------|
| Per-sample calling | `gatk HaplotypeCaller` | `-ERC GVCF`, `--native-pair-hmm-threads 4` | Per-sample gVCF |
| Database import | `gatk GenomicsDBImport` | `--batch-size 50`, `--reader-threads 4` | GenomicsDB workspace (temporary) |
| Joint genotyping | `gatk GenotypeGVCFs` | — | `cohort.vcf.gz` + index |

**WES:** Set `haplotypecaller.intervals` to capture kit BED.
**Joint calling:** `genomicsdb.intervals` is REQUIRED (`.list` format with genomic intervals).

**Script:** `scripts/run_variant_calling.sh` — Supports `--mode single` (per-sample gVCF) and `--mode joint` (cohort genotyping). Automatically indexes output VCFs.

---

### Stage 4b: Gene Quantification (RNA only)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-t` | `exon` | Count at exon level |
| `-g` | `gene_id` | Summarize to gene level |
| `--primary` | enabled | Only count primary alignments |
| `-s` | `0` | Strandedness: 0=unstranded, 1=stranded, 2=reverse |
| `-p --countReadPairs` | (PE only) | Count fragments, not individual reads |
| `-T` | 16 | Threads |

**Per-sample outputs** (temporary — deleted after merge):
- `{sample}.featureCounts.txt` — Raw gene counts
- `{sample}.featureCounts.txt.summary` — Assignment statistics

**Final output:** `gene_counts_matrix.tsv` — Gene × sample count matrix (genes in rows, samples in columns). Ready for DESeq2, edgeR, or limma-voom.

**Scripts:**
- `scripts/run_featurecounts.sh` — Supports `--mode single` and `--mode merge`. Auto-detects PE from `--paired` flag.
- `scripts/merge_counts.py` (inline in Snakemake) — Merges per-sample counts on `Geneid` column with outer join, fills missing values with 0.

---

### Stage 5: Ancestry Estimation

#### BAM Staging

Before ancestry analysis, all final BAMs are symlinked to a staging directory with sanitized names:
- Pipeline suffixes (`.sorted.markdup.recal`, etc.) are stripped
- Underscores are converted to hyphens (ADMIXTURE requirement)
- Only one BAM per sample is staged (avoids the `.sorted` suffix duplication issue)

#### Ancestry Container

| Component | Version | Role |
|-----------|---------|------|
| GATK | 3.8 | Genotype likelihood at ~10,000 AIM sites |
| PLINK | 1.9 | Merge with 1000 Genomes reference panel |
| ADMIXTURE | 1.3.0 | Supervised ancestry estimation |
| 1000 Genomes Phase 3 | 2,504 individuals | Reference panel (AFR, AMR, EAS, EUR, SAS) |

**Container execution:**
```bash
singularity run \
  --bind {staging_dir}:/input_data:ro \
  --bind {results_dir}/ancestry:/output \
  {container} \
  -bam_dir /input_data -outdir /output \
  -seq_type {rna|exome|wgs} -threads {threads} {--subpops}
```

**Seq-type mapping:** `rnaseq → rna`, `wes → exome`, `wgs → wgs`. Mixed cohorts default to `exome`.

#### Ancestry Visualization

**Script:** `scripts/plot_ancestry.R`

**Outputs:**
- PCA projection (4-panel: PC1-PC2, PC2-PC3, PC3-PC4, PC4-PC5)
- Combined ancestry barplot (K=5 super-populations + K=23 sub-populations)
- Individual barplots (K=5 only, K=23 only)
- All outputs in PNG (300 dpi) and SVG formats

**Color scheme:**
- AFR: `#E41A1C` (red), AMR: `#FF7F00` (orange), EAS: `#4DAF4A` (green), EUR: `#377EB8` (blue), SAS: `#984EA3` (purple)
- Sub-populations use gradient shades within each super-population

**Sample ordering:** Decreasing African (AFR) ancestry proportion.

**Adaptive figure sizing:**
- ≤15 samples: 0.7 in per sample + 4 in
- 15–50 samples: 0.45 in per sample + 4 in
- >50 samples: 0.3 in per sample + 6 in (capped at 48 in)

**Script:** `scripts/run_ancestry.sh` — Handles BAM staging, index detection (`.bai` or `.bam.bai`), container execution, and result verification.

---

### Stage 6: Reports

Three PDF reports generated using fpdf2 with optional lab logo:

| Report | Script | Input | Flags |
|--------|--------|-------|-------|
| Read QC | `scripts/fastq_qc_report.py` | fastp JSON | Pass rate < 80%, Q30 < 80%, Dup > 30% |
| Alignment QC | `scripts/alignment_qc_report.py` | flagstat + dup_metrics | Mapping < 90%, Dup > 30% |
| Ancestry | `scripts/ancestry_qc_report.py` | ancestry_summary TSV + barplot PNG | Dominant ancestry highlighted |

Each report produces both a machine-readable TSV summary and a visual PDF.

---

## Execution Profiles

### SLURM Profile (`profiles/slurm/config.yaml`)

```yaml
executor: slurm
default-resources:
  slurm_account: "yateslab_genomics"
  slurm_partition: "cpu"      # cluster-specific — check `sinfo -s`
  mem_mb: 8000
  runtime: 60
  cpus_per_task: 1
  slurm_extra: "'--mail-type=FAIL'"
jobs: 50                    # Max concurrent SLURM jobs
use-envmodules: true
latency-wait: 120           # Seconds to wait for NFS file propagation
keep-going: true
rerun-incomplete: true
printshellcmds: true
```

**Per-rule resource overrides (key rules):**

| Rule | Threads | Memory | Runtime |
|------|---------|--------|---------|
| fastp_pe / fastp_se | 4 | 8 GB | 120 min |
| bwa_mem2_align | 16 | 48 GB | 480 min |
| hisat2_align | 16 | 32 GB | 480 min |
| samtools_sort | 8 | 32 GB | 180 min |
| mark_duplicates | 1 | 32 GB | 240 min |
| base_recalibrator / apply_bqsr | 1 | 16 GB | 360 min |
| haplotype_caller | 4 | 16 GB | 720 min |
| genomics_db_import | 4 | 48 GB | 480 min |
| genotype_gvcfs | 4 | 32 GB | 480 min |
| run_ancestry | 4 | 80 GB | 1440 min |

**Important:** Run Snakemake from the **login node** when using the SLURM profile. Running from within a SLURM job (compute node) can cause job submission failures.

### Local Profile (`profiles/local/config.yaml`)

```yaml
cores: 30
use-envmodules: true
latency-wait: 30
keep-going: true
rerun-incomplete: true
resources:
  mem_mb: 64000
```

Suitable for interactive sessions on a single compute node.

### Running the Pipeline

```bash
# Dry run (check DAG without executing)
snakemake --profile profiles/slurm -n

# Full execution via SLURM (from login node)
snakemake --profile profiles/slurm

# Full execution locally (from compute node)
snakemake --profile profiles/local

# Resume after failure
snakemake --profile profiles/slurm --rerun-incomplete
snakemake --profile profiles/local --rerun-incomplete
```

### Re-running Specific Stages

Use `--forcerun` to re-execute specific rules without re-running the entire pipeline. Snakemake will also re-run any downstream rules that depend on the forced rules. Always dry-run first with `-n`.

#### Examples

```bash
# Re-run a single sample's alignment
snakemake --profile profiles/local \
  --forcerun hisat2_align \
  results/alignment/bam/SampleA.raw.bam

# Re-run featureCounts + merge for all samples
snakemake --profile profiles/local \
  --forcerun featurecounts merge_counts

# Re-run ancestry (container + plots + report)
snakemake --profile profiles/local \
  --forcerun run_ancestry plot_ancestry ancestry_qc_report

# Re-run ancestry plots only (after updating plot_ancestry.R)
snakemake --profile profiles/local \
  --forcerun plot_ancestry ancestry_qc_report

# Re-run full ancestry from staging (full re-estimation)
snakemake --profile profiles/local \
  --forcerun stage_bams_for_ancestry run_ancestry plot_ancestry ancestry_qc_report

# Re-run all reports only
snakemake --profile profiles/local \
  --forcerun fastq_qc_report alignment_qc_report ancestry_qc_report

# Re-run MultiQC only
snakemake --profile profiles/local \
  --forcerun multiqc

# Reduce threads for a specific rule (e.g., large sample hitting memory)
snakemake --profile profiles/local \
  --set-threads hisat2_align=8 \
  --rerun-incomplete

# Dry-run first — preview what will be re-run before executing
snakemake --profile profiles/local \
  --forcerun featurecounts merge_counts -n
```

#### Available Rule Names for `--forcerun`

| Stage | Rule Names |
|-------|-----------|
| **QC** | `fastp_pe`, `fastp_se`, `fastqc_raw_r1`, `fastqc_raw_r2`, `fastqc_trimmed_r1`, `fastqc_trimmed_r2`, `multiqc` |
| **Alignment** | `bwa_mem2_align` (DNA), `hisat2_align` (RNA), `samtools_sort`, `mark_duplicates`, `samtools_index`, `samtools_flagstat` |
| **BQSR** | `base_recalibrator`, `apply_bqsr`, `index_recal_bam` |
| **Variant Calling** | `haplotype_caller`, `genomics_db_import`, `genotype_gvcfs` |
| **Quantification** | `featurecounts`, `merge_counts` |
| **Ancestry** | `stage_bams_for_ancestry`, `run_ancestry`, `plot_ancestry` |
| **Reports** | `fastq_qc_report`, `alignment_qc_report`, `ancestry_qc_report` |

---

## Output Directory Structure

Files marked 🗑️ are **temporary** — Snakemake auto-deletes them after downstream rules consume them.

```
results/
├── qc/
│   ├── fastp/
│   │   ├── {sample}_R1.trimmed.fastq.gz     🗑️  Deleted after alignment
│   │   ├── {sample}_R2.trimmed.fastq.gz     🗑️  Deleted after alignment (PE)
│   │   ├── {sample}_fastp.json               Kept (MultiQC + reports)
│   │   └── {sample}_fastp.html              🗑️  Redundant with MultiQC
│   ├── fastqc_raw/
│   │   ├── {sample}_{R1,R2}_fastqc.html      Kept
│   │   └── {sample}_{R1,R2}_fastqc.zip      🗑️  Deleted after MultiQC
│   ├── fastqc_trimmed/
│   │   ├── {sample}_{R1,R2}.trimmed_fastqc.html  Kept
│   │   └── {sample}_{R1,R2}.trimmed_fastqc.zip  🗑️
│   ├── multiqc/multiqc_report.html           Aggregated QC report
│   └── reports/
│       ├── fastq_qc_summary.tsv
│       └── fastq_qc_report.pdf
│
├── alignment/
│   ├── bam/
│   │   ├── {sample}.raw.bam                 🗑️  Deleted after sorting
│   │   ├── {sample}.sorted.bam              🗑️  Deleted after MarkDuplicates
│   │   ├── {sample}.sorted.markdup.bam       RNA final BAM
│   │   ├── {sample}.sorted.markdup.bam.bai
│   │   ├── {sample}.sorted.markdup.recal.bam DNA final BAM (if BQSR)
│   │   └── {sample}.sorted.markdup.recal.bam.bai
│   ├── bqsr/{sample}.recal_data.table        BQSR model (DNA only)
│   ├── metrics/
│   │   ├── {sample}.dup_metrics.txt
│   │   ├── {sample}.flagstat.txt
│   │   └── {sample}.hisat2_summary.txt       RNA only
│   └── reports/
│       ├── alignment_qc_summary.tsv
│       └── alignment_qc_report.pdf
│
├── variants/                                  DNA only
│   ├── gvcf/{sample}.g.vcf.gz (+.tbi)        Per-sample gVCFs (kept)
│   ├── genomicsdb/                           🗑️  Deleted after genotyping
│   └── joint/
│       ├── cohort.vcf.gz                      Final joint VCF
│       └── cohort.vcf.gz.tbi
│
├── counts/                                    RNA only
│   ├── {sample}.featureCounts.txt            🗑️  Deleted after merge
│   ├── {sample}.featureCounts.txt.summary    🗑️
│   └── gene_counts_matrix.tsv                Final count matrix
│
├── ancestry/
│   ├── .staging/                             🗑️  Staging symlinks
│   ├── ancestry_summary_superpops.tsv         K=5 continental ancestry
│   ├── ancestry_summary_subpops.tsv           K=23 (if subpops enabled)
│   ├── ancestry_pca_coordinates.tsv           PCA projection
│   ├── ancestry_barplot_combined.png/svg      Combined barplot
│   ├── ancestry_barplot_superpops.png/svg     K=5 barplot
│   ├── ancestry_barplot_subpops.png/svg       K=23 barplot (if enabled)
│   ├── {sample}.ancestry_superpops.tsv        Per-sample ancestry
│   └── reports/ancestry_report.pdf
│
└── logs/                                      Per-rule, per-sample logs
    ├── fastp/{sample}.log
    ├── fastqc/{sample}_{R1,R2}_{raw,trimmed}.log
    ├── bwa_mem2/{sample}.log
    ├── hisat2/{sample}.log
    ├── samtools_sort/{sample}.log
    ├── markdup/{sample}.log
    ├── bqsr/{sample}.{base_recalibrator,apply_bqsr}.log
    ├── haplotypecaller/{sample}.log
    ├── featurecounts/{sample}.log
    ├── ancestry/ancestry_batch.log
    └── reports/
```

---

## Storage Optimization

The pipeline uses Snakemake's `temp()` directive to auto-delete intermediate files, saving **10–20 TB per 50-sample cohort**:

| Intermediate File | Deleted After | Savings Per Sample |
|-------------------|---------------|-------------------|
| Trimmed FASTQs | Alignment completes | 50–200 GB |
| Raw BAMs | Sorting completes | 50–200 GB |
| Sorted BAMs | MarkDuplicates completes | 50–200 GB |
| fastp HTML | MultiQC completes | 2–8 MB |
| FastQC ZIP archives | MultiQC completes | ~500 KB |
| GenomicsDB workspace | Joint genotyping completes | 50–200 GB total |
| Per-sample featureCounts | Count merge completes | 1–10 MB |
| Ancestry staging dir | Ancestry run completes | Symlinks only |

**What is kept permanently:** fastp JSON, FastQC HTML, MultiQC report, final BAMs, gVCFs, joint VCF, merged count matrix, ancestry results, all PDF reports, all logs.

---

## Standalone Script Usage

Every pipeline step has a standalone Bash script that can be run independently:

| Script | Usage |
|--------|-------|
| `generate_samples.py` | `python scripts/generate_samples.py --input-dir /data/ --seq-type rnaseq --output config/samples.tsv` |
| `run_fastp.sh` | `bash scripts/run_fastp.sh --r1 R1.fq.gz [--r2 R2.fq.gz] --out-r1 ... --json ... --html ... --qual 20 --minlen 50 --threads 16` |
| `run_fastqc.sh` | `bash scripts/run_fastqc.sh --input R1.fq.gz --outdir qc/ --threads 8` |
| `run_multiqc.sh` | `bash scripts/run_multiqc.sh --indir qc/ --outdir multiqc/` |
| `run_bwa_mem2.sh` | `bash scripts/run_bwa_mem2.sh --r1 R1.fq.gz --r2 R2.fq.gz --ref ref.fa --sample S1 --outdir align/ --seq-type WES --threads 16` |
| `run_hisat2.sh` | `bash scripts/run_hisat2.sh --r1 R1.fq.gz [--r2 R2.fq.gz] --index grch38 --sample S1 --outdir align/ --threads 16` |
| `run_bqsr.sh` | `bash scripts/run_bqsr.sh --bam S1.markdup.bam --ref ref.fa --sample S1 --outdir align/ --known-sites dbsnp.vcf.gz --known-sites mills.vcf.gz` |
| `run_variant_calling.sh` | `bash scripts/run_variant_calling.sh --mode single --bam S1.recal.bam --ref ref.fa --sample S1 --outdir variants/` |
| `run_featurecounts.sh` | `bash scripts/run_featurecounts.sh --mode single --bam S1.markdup.bam --gtf anno.gtf --sample S1 --outdir counts/ [--paired] --strandedness 0` |
| `run_ancestry.sh` | `bash scripts/run_ancestry.sh --bam-dir bams/ --outdir ancestry/ --container pipeline.sif --seq-type rna [--subpops] --threads 30` |
| `plot_ancestry.R` | `Rscript scripts/plot_ancestry.R {results_dir} [prefix] [metadata.tsv]` |

All Bash scripts accept `--help` for full usage documentation.

---

## Resource Requirements

| Stage | RAM | Threads | Time (per sample) | Disk (per sample) |
|-------|-----|---------|-------------------|-------------------|
| fastp | 8 GB | 4–16 | 5–15 min | ~same as input |
| FastQC | 4 GB | 2–8 | 5–10 min | ~2 MB |
| bwa-mem2 alignment | 48 GB | 16 | 30–90 min | ~2× input size |
| HISAT2 alignment | 32 GB | 16 | 20–60 min | ~1.5× input size |
| samtools sort | 32 GB | 8 | 15–30 min | ~same as BAM |
| MarkDuplicates | 32 GB | 1 | 30–60 min | ~same as BAM |
| BQSR (2-step) | 16 GB | 1 | 60–120 min | ~same as BAM |
| HaplotypeCaller | 16 GB | 4 | 2–12 hours | ~0.5–5 GB (gVCF) |
| GenomicsDBImport | 48 GB | 4 | 1–4 hours | ~50–200 GB |
| GenotypeGVCFs | 32 GB | 4 | 1–4 hours | ~0.1–1 GB (VCF) |
| featureCounts | 8 GB | 4–16 | 5–15 min | ~10 MB |
| Ancestry | 80 GB | 4–30 | 15–60 min | ~50 MB |

---

## Troubleshooting

### Common Issues

| Error | Cause | Fix |
|-------|-------|-----|
| `Cannot load module "multiqc/1.20"` | MultiQC requires python311 pre-loaded | Pipeline handles this; ensure `python311` module exists |
| `Pillow not available - fpdf2 cannot insert images` | Missing Python package | `pip install --user Pillow` |
| `Invalid SLURM account` | Wrong account name in profile | Update `slurm_account` in `profiles/slurm/config.yaml` |
| `Directory cannot be locked` | Previous Snakemake run crashed | `snakemake --profile profiles/slurm --unlock` |
| `executor: invalid choice: 'slurm'` | Missing SLURM plugin | `pip install snakemake-executor-plugin-slurm` |
| SLURM jobs fail from compute node | Submitting SLURM jobs within a SLURM job | Run Snakemake from login node, or use local profile |
| Duplicate ancestry results | BAM filename not properly sanitized | Fixed: staging now strips all pipeline suffixes |
| `--countReadPairs` not recognized | Subread version < 2.0 | Update to Subread >= 2.0.1 |

### Failure Recovery

The pipeline failure banner shows the correct re-run command:
```
============================================================
  PIPELINE FAILED — check logs in results/logs/
  Re-run with:
    SLURM:  snakemake --profile profiles/slurm --rerun-incomplete
    Local:  snakemake --profile profiles/local --rerun-incomplete
============================================================
```

---

## Authors

Samuel Mwamburi — Yates Lab, Johns Hopkins University
