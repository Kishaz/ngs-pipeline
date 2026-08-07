# NGS Pipeline

A production-grade Snakemake pipeline for processing **RNAseq**, **WES**, and **WGS** data through quality control, alignment, base quality recalibration, variant calling / gene quantification, and genetic ancestry estimation.

Supports both **paired-end** and **single-end** reads, **FASTQ** and **pre-aligned BAM** inputs with automatic sample discovery. Each stage writes outputs to dedicated directories in standard formats for downstream analysis.

---

## Pipeline Overview

```
Raw FASTQ
  │
  ▼
┌──────────────────────────────────────────────────────┐
│  STAGE 1: READ QC                                    │
│  fastp → FastQC (raw + trimmed) → MultiQC            │
│  Output: trimmed FASTQs, QC reports, PDF summary     │
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

## Setup

A condensed on-ramp — see [Quick Start](#quick-start) for full detail and [Tools and Versions Used](#tools-and-versions-used) for the exact environment this pipeline was run with.

1. **Dependencies** — Snakemake (≥ 8.0, typically in a conda env), plus the report/plugin packages:
   ```bash
   pip install --user fpdf2 pandas Pillow
   pip install snakemake-executor-plugin-slurm   # SLURM profile only
   ```
2. **References** (one-time, ~60 GB) — `setup_references.sh` downloads and indexes GRCh38, GENCODE v48, the GATK known-sites bundle, and the ancestry container:
   ```bash
   bash scripts/setup_references.sh --outdir /path/to/ngs_resources
   ```
3. **Configure** — edit `config/config.yaml`: resource paths, `input_dir` (or `config/samples.tsv`), and your cluster's `module`/`conda` tool names.
4. **Run**:
   ```bash
   snakemake --profile profiles/slurm      # cluster;  profiles/local for an interactive session
   ```

---

## Pipeline Stages — Detailed

### Stage 1: Read Quality Control

| Item | Detail |
|------|--------|
| **Tool** | [fastp](https://github.com/OpenGene/fastp) |
| **Input** | Raw FASTQ files (`.fastq.gz`), PE or SE |
| **Output** | Trimmed FASTQs *(temporary — auto-deleted after alignment)*, per-sample JSON reports |
| **Parameters** | `--qualified_quality_phred 20`, `--length_required 50`, `--detect_adapter_for_pe` (PE only) |
| **Role** | Removes adapters, trims low-quality bases, filters short reads. Ensures only high-quality reads enter alignment. |
| **Script** | `scripts/run_fastp.sh` |

| Item | Detail |
|------|--------|
| **Tool** | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| **Input** | Raw FASTQs (pre-trim) and trimmed FASTQs (post-trim) |
| **Output** | Per-file HTML quality reports *(ZIP archives auto-deleted after use)* |
| **Role** | Assesses per-base quality, GC content, adapter contamination, duplication levels, and overrepresented sequences. Run before and after trimming to verify improvement. |
| **Script** | `scripts/run_fastqc.sh` |

| Item | Detail |
|------|--------|
| **Tool** | [MultiQC](https://multiqc.info/) |
| **Input** | All FastQC and fastp reports |
| **Output** | Single aggregated HTML report |
| **Role** | Combines QC metrics from all samples into one interactive report for batch-level assessment. |
| **Script** | `scripts/run_multiqc.sh` |

| Item | Detail |
|------|--------|
| **Tool** | Custom Python (fpdf2) |
| **Input** | fastp JSON reports |
| **Output** | `fastq_qc_summary.tsv` + `fastq_qc_report.pdf` |
| **Role** | Generates a shareable PDF with per-sample metrics table. Flags samples with pass rate < 80%, Q30 < 80%, or duplication > 30%. |
| **Script** | `scripts/fastq_qc_report.py` |

---

### Stage 2: Sequence Alignment

**DNA samples (WES/WGS):**

| Item | Detail |
|------|--------|
| **Tool** | [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) |
| **Input** | Trimmed FASTQs + GRCh38 reference (bwa-mem2 indexed) |
| **Output** | Unsorted BAM (piped to samtools) |
| **Parameters** | `-Y` (soft-clip supplementary), `-t 16`, full `@RG` header |
| **Role** | SIMD-accelerated short-read aligner. Maps reads to the reference genome with read group metadata for downstream GATK compatibility. |
| **Script** | `scripts/run_bwa_mem2.sh` |

**RNA samples:**

| Item | Detail |
|------|--------|
| **Tool** | [HISAT2](http://daehwankimlab.github.io/hisat2/) |
| **Input** | Trimmed FASTQs + HISAT2 splice-aware index |
| **Output** | Unsorted BAM + alignment summary |
| **Parameters** | `--dta` (downstream transcript assembly), `-p 16`, full `--rg` header |
| **Role** | Splice-aware aligner for RNA-seq. Uses a graph-based index with known splice sites from GENCODE annotation for accurate junction mapping. |
| **Script** | `scripts/run_hisat2.sh` |

**Shared post-alignment (all samples):**

| Step | Tool | Input | Output | Role |
|------|------|-------|--------|------|
| Sort | `samtools sort` | Raw BAM | Coordinate-sorted BAM | Orders reads by genomic position for efficient downstream access. |
| MarkDuplicates | `gatk MarkDuplicates` | Sorted BAM | Deduplicated BAM + metrics | Flags PCR/optical duplicates to prevent variant calling bias. Duplicates are flagged, not removed. |
| Index | `samtools index` | Markdup BAM | `.bai` index | Enables random access to BAM regions. Required by GATK, IGV, and most downstream tools. |
| Flagstat | `samtools flagstat` | Markdup BAM | `flagstat.txt` | Computes alignment summary (total reads, mapped %, properly paired %, etc.). |

---

### Stage 3: Base Quality Score Recalibration (BQSR)

> **DNA samples only.** RNA samples skip this stage entirely.

| Item | Detail |
|------|--------|
| **Tool** | [GATK BaseRecalibrator + ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR) |
| **Input** | Markdup BAM + reference FASTA + known variant sites (dbSNP, Mills indels, known indels) |
| **Output** | Recalibration table (`.recal_data.table`) + recalibrated BAM (`.sorted.markdup.recal.bam`) |
| **Role** | Corrects systematic errors in base quality scores assigned by the sequencer. Machine learning model learns error patterns from positions matching known variants, then adjusts all base qualities. Produces the analysis-ready BAM used for variant calling. |
| **Script** | `scripts/run_bqsr.sh` |

**Known sites files (from GATK resource bundle):**
- `Homo_sapiens_assembly38.dbsnp138.vcf.gz` — common human variants
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` — validated insertion/deletion sites
- `Homo_sapiens_assembly38.known_indels.vcf.gz` — additional known indels

**Why BQSR matters:** Without recalibration, base quality scores may systematically overestimate or underestimate confidence, leading to false positive or false negative variant calls. BQSR is a standard step in GATK best practices for germline variant discovery.

---

### Stage 4a: Variant Calling (DNA only)

| Step | Tool | Input | Output | Role |
|------|------|-------|--------|------|
| Per-sample calling | `gatk HaplotypeCaller` | Recalibrated BAM + reference | Per-sample gVCF (`.g.vcf.gz`) | Local de novo assembly of haplotypes. Calls SNPs and indels with likelihoods for all possible genotypes. `-ERC GVCF` mode emits reference confidence blocks for joint genotyping. |
| Database import | `gatk GenomicsDBImport` | All sample gVCFs + intervals | GenomicsDB workspace | Consolidates per-sample gVCFs into an efficient columnar database for joint analysis. Scales to thousands of samples. |
| Joint genotyping | `gatk GenotypeGVCFs` | GenomicsDB + reference | Joint cohort VCF (`cohort.vcf.gz`) | Re-genotypes all samples simultaneously, leveraging population-level information to improve accuracy at rare variant sites. Produces the final multi-sample VCF. |

**Script:** `scripts/run_variant_calling.sh`

**WES note:** Set `variant_calling.haplotypecaller.intervals` to your capture kit BED file to restrict calling to on-target regions (faster, fewer artifacts).

---

### Stage 4b: Gene Quantification (RNA only)

| Item | Detail |
|------|--------|
| **Tool** | [featureCounts](http://subread.sourceforge.net/) (from Subread) |
| **Input** | Markdup BAM + GENCODE GTF annotation |
| **Output** | Merged gene counts matrix (`gene_counts_matrix.tsv`) *(per-sample counts are temporary — auto-deleted after merging)* |
| **Parameters** | `-t exon -g gene_id --primary`, strandedness from config |
| **Role** | Assigns aligned reads to genomic features (exons → genes). The merged matrix has genes as rows and samples as columns — ready for DESeq2, edgeR, or limma-voom. |
| **Script** | `scripts/run_featurecounts.sh` |

---

### Stage 5: Ancestry Estimation

| Item | Detail |
|------|--------|
| **Tool** | Containerized pipeline: GATK 3.8 + PLINK 1.9 + ADMIXTURE 1.3.0 |
| **Input** | Final BAMs (all samples, symlinked to staging directory — staging auto-deleted after ancestry completes) |
| **Output** | Per-sample ancestry proportions, combined summary TSV, PCA projection, barplots (PNG + SVG) |
| **Role** | Estimates continental ancestry composition for each sample using supervised ADMIXTURE with the 1000 Genomes Phase 3 reference panel (2,283 individuals, 5 super-populations). |
| **Script** | `scripts/run_ancestry.sh` (staging + Singularity) + `scripts/plot_ancestry.R` (visualization with optional metadata annotations) |

**Processing steps inside the container:**
1. GATK 3.8 pileup at ~10,000 ancestry-informative marker sites
2. Likelihood-based genotype calling from pileup data
3. PLINK merge with 1000 Genomes reference panel
4. ADMIXTURE supervised estimation (K=5 continental, optional K=23 sub-population)
5. PCA projection onto reference coordinate space

**Output super-populations:** AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian)

---

### Stage 6: Reports

Three PDF reports are generated, one per major stage:

| Report | Input | Content |
|--------|-------|---------|
| **Read QC Report** | fastp JSON files | Per-sample: raw reads, clean reads, pass rate, Q30, GC%, duplication. Flags low-quality samples in red. |
| **Alignment QC Report** | flagstat + MarkDuplicates metrics | Per-sample: total reads, mapped reads, mapping rate, properly paired %, duplication rate. Flags mapping < 90% or dup > 30%. |
| **Ancestry Report** | Superpop TSVs + combined barplot | Page 1: ancestry proportions table with dominant population highlighted. Page 2: stacked barplot figure. |

---

## Output Directory Structure

Files marked 🗑️ are **temporary** — Snakemake auto-deletes them once all downstream rules have consumed them. This saves **10-20 TB per 50-sample cohort**.

```
results/
├── qc/
│   ├── fastqc_raw/                     FastQC HTML reports on raw reads
│   ├── fastp/
│   │   ├── {sample}_fastp.json          QC metrics (kept — used by MultiQC + reports)
│   │   ├── {sample}_fastp.html          🗑️ Per-sample HTML (redundant with MultiQC)
│   │   ├── {sample}_R1.trimmed.fq.gz   🗑️ Trimmed reads (deleted after alignment)
│   │   └── {sample}_R2.trimmed.fq.gz   🗑️
│   ├── fastqc_trimmed/                 FastQC HTML reports on trimmed reads
│   ├── multiqc/multiqc_report.html     Aggregated QC report
│   └── reports/
│       ├── fastq_qc_summary.tsv        Machine-readable QC metrics
│       └── fastq_qc_report.pdf         Shareable PDF report
│
├── alignment/
│   ├── bam/
│   │   ├── {sample}.sorted.markdup.bam          Deduplicated BAM (RNA final BAM)
│   │   ├── {sample}.sorted.markdup.bam.bai      BAM index
│   │   ├── {sample}.sorted.markdup.recal.bam    Recalibrated BAM (DNA final BAM)
│   │   └── {sample}.sorted.markdup.recal.bam.bai
│   ├── bqsr/
│   │   └── {sample}.recal_data.table             BQSR recalibration model
│   ├── metrics/
│   │   ├── {sample}.dup_metrics.txt              MarkDuplicates metrics
│   │   ├── {sample}.flagstat.txt                 Alignment statistics
│   │   └── {sample}.hisat2_summary.txt           HISAT2 summary (RNA only)
│   └── reports/
│       ├── alignment_qc_summary.tsv
│       └── alignment_qc_report.pdf
│
├── variants/                            (DNA only)
│   ├── gvcf/{sample}.g.vcf.gz          Per-sample gVCFs (kept for re-genotyping)
│   ├── genomicsdb/                      🗑️ GenomicsDB workspace (deleted after genotyping)
│   └── joint/
│       ├── cohort.vcf.gz               Joint-called VCF
│       └── cohort.vcf.gz.tbi           VCF index
│
├── counts/                              (RNA only)
│   ├── {sample}.featureCounts.txt       🗑️ Per-sample counts (deleted after merge)
│   └── gene_counts_matrix.tsv          Combined counts matrix (final deliverable)
│
├── ancestry/
│   ├── .staging/                        🗑️ Symlinked BAMs (deleted after ancestry)
│   ├── {sample}.ancestry_superpops.tsv  Per-sample continental ancestry
│   ├── ancestry_summary_superpops.tsv   Combined summary
│   ├── ancestry_pca.pdf                 PCA projection plot
│   ├── ancestry_barplot_combined.png    Stacked ancestry barplot
│   └── reports/ancestry_report.pdf      Shareable PDF report
│
└── logs/                                Per-rule, per-sample log files
    ├── fastp/
    ├── fastqc/
    ├── bwa_mem2/
    ├── hisat2/
    ├── samtools_sort/
    ├── markdup/
    ├── bqsr/
    ├── haplotypecaller/
    ├── featurecounts/
    ├── ancestry/
    └── reports/
```

### Intermediate files auto-deleted by Snakemake

| File | Deleted after | Size saved |
|------|---------------|------------|
| Trimmed FASTQs | Alignment completes | 50-200 GB/sample |
| Raw BAMs | Sorting completes | 50-200 GB/sample |
| Sorted BAMs | MarkDuplicates completes | 50-200 GB/sample |
| fastp HTML | MultiQC completes | 2-8 MB/sample |
| FastQC ZIP archives | MultiQC completes | ~500 KB/sample |
| GenomicsDB workspace | Joint genotyping completes | 50-200 GB total |
| Per-sample featureCounts | Merge completes | 1-10 MB/sample |
| Ancestry staging dir | Ancestry run completes | Symlinks only |

---

## Quick Start

### 1. Install Python dependencies and SLURM plugin

```bash
pip install --user fpdf2 pandas Pillow
pip install snakemake-executor-plugin-slurm   # required for SLURM profile
```

### 2. Set up resources

```bash
# On HPC: download and index all references (~60 GB, ~4 hours)
sbatch --mem=200G --cpus-per-task=16 --time=06:00:00 \
    scripts/setup_references.sh --outdir /path/to/ngs_resources
```

### 3. Configure input

Edit `config/config.yaml`:

```yaml
# Option A: Auto-discover samples from a directory (recommended)
input_dir: "/path/to/your/fastqs_or_bams"
default_seq_type: "rnaseq"    # rnaseq | wes | wgs

# Option B: Manual sample manifest (leave input_dir empty)
input_dir: ""
samples: "config/samples.tsv"
```

**Option A — Auto-discovery** scans `input_dir` for `.fastq.gz` and `.bam` files on every run:
- Detects paired-end vs single-end from standard Illumina naming (`_R1/_R2`, `_1/_2`)
- Auto-detects BAM processing stage from filename (`.sorted.markdup.bam`, etc.)
- Generates `config/samples.tsv` automatically

**Option B — Manual manifest** uses a hand-crafted `config/samples.tsv` (tab-delimited):

```
sample      R1                              R2                              seq_type
PatientA    /data/PatientA_R1.fastq.gz      /data/PatientA_R2.fastq.gz      wes
PatientB    /data/PatientB_R1.fastq.gz      /data/PatientB_R2.fastq.gz      wgs
PatientC    /data/PatientC_R1.fastq.gz                                      rnaseq
```

- `R2` is empty for single-end reads
- `seq_type`: `wes`, `wgs`, or `rnaseq`

### 4. BAM Input Mode

Pre-aligned BAMs can enter the pipeline at the appropriate stage. The processing stage is auto-detected from the filename:

| Filename pattern | Detected stage | Pipeline starts at |
|-----------------|----------------|-------------------|
| `*.sorted.markdup.recal.bam` | `recal` | Indexing only |
| `*.sorted.markdup.bam` | `markdup` | Indexing + flagstat |
| `*.sorted.bam` | `sorted` | MarkDuplicates |
| `*.raw.bam` or `*.bam` | `raw` | samtools sort |

BAM samples skip QC (fastp/FastQC). For manual TSV, add `bam` and `bam_stage` columns:

```
sample      R1      R2      seq_type    bam                                     bam_stage
PatientE                    wes         /data/PatientE.sorted.markdup.bam       markdup
```

### 5. Update config

Edit `config/config.yaml` — set paths to your resources, module names for your cluster, and project name.

### 6. Dry-run

```bash
snakemake --profile profiles/slurm -n
```

### 7. Execute

```bash
# SLURM (job submission)
snakemake --profile profiles/slurm

# Local (interactive session)
snakemake --profile profiles/local
```

### 8. Resume after failure

The pipeline prints the correct re-run command on failure:

```bash
# SLURM
snakemake --profile profiles/slurm --rerun-incomplete

# Local
snakemake --profile profiles/local --rerun-incomplete
```

### 9. Re-running specific stages

Use `--forcerun` to re-execute specific rules without re-running the entire pipeline. Snakemake will also re-run any downstream rules that depend on the forced rules.

**Re-run a single sample's alignment:**
```bash
snakemake --profile profiles/local \
  --forcerun hisat2_align \
  results/alignment/bam/SampleA.raw.bam
```

**Re-run featureCounts + merge for all samples:**
```bash
snakemake --profile profiles/local \
  --forcerun featurecounts merge_counts
```

**Re-run ancestry (container + plots + report):**
```bash
snakemake --profile profiles/local \
  --forcerun run_ancestry plot_ancestry ancestry_qc_report
```

**Re-run ancestry plots only (after updating `plot_ancestry.R`):**
```bash
snakemake --profile profiles/local \
  --forcerun plot_ancestry ancestry_qc_report
```

**Re-run all reports only:**
```bash
snakemake --profile profiles/local \
  --forcerun fastq_qc_report alignment_qc_report ancestry_qc_report
```

**Re-run from staging through ancestry (full ancestry re-estimation):**
```bash
snakemake --profile profiles/local \
  --forcerun stage_bams_for_ancestry run_ancestry plot_ancestry ancestry_qc_report
```

**Re-run MultiQC only:**
```bash
snakemake --profile profiles/local \
  --forcerun multiqc
```

**Dry-run first** — always add `-n` to preview what will be re-run before executing:
```bash
snakemake --profile profiles/local \
  --forcerun featurecounts merge_counts -n
```

**Reduce threads for a specific rule** (e.g., large sample hitting memory limits):
```bash
snakemake --profile profiles/local \
  --set-threads hisat2_align=8 \
  --rerun-incomplete
```

### Available rule names for `--forcerun`

| Stage | Rule names |
|-------|-----------|
| QC | `fastp_pe`, `fastp_se`, `fastqc_raw_r1`, `fastqc_raw_r2`, `fastqc_trimmed_r1`, `fastqc_trimmed_r2`, `multiqc` |
| Alignment | `bwa_mem2_align` (DNA), `hisat2_align` (RNA), `samtools_sort`, `mark_duplicates`, `samtools_index`, `samtools_flagstat` |
| BQSR | `base_recalibrator`, `apply_bqsr`, `index_recal_bam` |
| Variant calling | `haplotype_caller`, `genomics_db_import`, `genotype_gvcfs` |
| Quantification | `featurecounts`, `merge_counts` |
| Ancestry | `stage_bams_for_ancestry`, `run_ancestry`, `plot_ancestry` |
| Reports | `fastq_qc_report`, `alignment_qc_report`, `ancestry_qc_report` |

---

## Standalone Script Usage

Every pipeline step has a standalone Bash script that can be run independently outside Snakemake:

| Script | Usage |
|--------|-------|
| `generate_samples.py` | `python scripts/generate_samples.py --input-dir /data/ --seq-type rnaseq --output config/samples.tsv` |
| `run_fastp.sh` | `bash scripts/run_fastp.sh --r1 R1.fq.gz [--r2 R2.fq.gz] --out-r1 ... --json ... --html ...` |
| `run_fastqc.sh` | `bash scripts/run_fastqc.sh --input R1.fq.gz --outdir qc/` |
| `run_multiqc.sh` | `bash scripts/run_multiqc.sh --indir qc/ --outdir multiqc/` |
| `run_bwa_mem2.sh` | `bash scripts/run_bwa_mem2.sh --r1 R1.fq.gz --r2 R2.fq.gz --ref ref.fa --sample S1 --outdir align/ --seq-type WES` |
| `run_hisat2.sh` | `bash scripts/run_hisat2.sh --r1 R1.fq.gz --r2 R2.fq.gz --index grch38 --sample S1 --outdir align/` |
| `run_bqsr.sh` | `bash scripts/run_bqsr.sh --bam S1.markdup.bam --ref ref.fa --sample S1 --outdir align/ --known-sites dbsnp.vcf.gz --known-sites mills.vcf.gz` |
| `run_variant_calling.sh` | `bash scripts/run_variant_calling.sh --mode single --bam S1.recal.bam --ref ref.fa --sample S1 --outdir variants/` |
| `run_featurecounts.sh` | `bash scripts/run_featurecounts.sh --mode single --bam S1.markdup.bam --gtf anno.gtf --sample S1 --outdir counts/ --paired` |
| `run_ancestry.sh` | `bash scripts/run_ancestry.sh --bam-dir bams/ --outdir ancestry/ --container pipeline.sif --subpops` |
| `plot_ancestry.R` | `Rscript scripts/plot_ancestry.R ./results cohort metadata.csv -c race,gender` |

All scripts accept `--help` for full usage.

---

## Ancestry Visualization (`plot_ancestry.R`)

Standalone companion script that generates publication-quality ancestry figures from pipeline output.

### Usage

```bash
Rscript scripts/plot_ancestry.R <results_dir> [prefix] [metadata_file] [-c col1,col2,...]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `results_dir` | Yes | Path to ancestry pipeline output directory |
| `prefix` | No | Output filename prefix (default: `ancestry`) |
| `metadata_file` | No | Sample metadata file (CSV, TSV, or XLSX — auto-detected) |
| `-c COLUMNS` | No | Comma-separated column names to show as annotation strips (default: all columns) |

### Examples

```bash
# Basic — no metadata annotations
Rscript scripts/plot_ancestry.R ./results/ancestry PDAC

# With metadata — show all columns as annotation strips
Rscript scripts/plot_ancestry.R ./results/ancestry PDAC metadata.csv

# Show only race annotation
Rscript scripts/plot_ancestry.R ./results/ancestry PDAC metadata.csv -c race

# Show race and treatment annotations
Rscript scripts/plot_ancestry.R ./results/ancestry PDAC metadata.xlsx -c race,treatment

# Show race, gender, and grade
Rscript scripts/plot_ancestry.R ./results/ancestry PDAC metadata.csv -c race,gender,grade
```

### Metadata file format

Any tabular file with a `sample_id` column. Format is auto-detected:

| Extension | Reader |
|-----------|--------|
| `.xlsx` / `.xls` | `readxl::read_excel()` (requires `readxl` package) |
| `.tsv` / tab-delimited | `read.delim()` |
| `.csv` / comma-delimited | `read.csv()` |

The `-c` flag is **case-insensitive** (`Race`, `race`, `RACE` all work). Invalid column names are warned and skipped.

### Output files

| File | Description |
|------|-------------|
| `{prefix}_pca.png/.svg` | PCA projection (4 panels: PC1-5 pairs) |
| `{prefix}_barplot_combined.png/.svg` | Combined figure: metadata strips + K=5 + K=23 barplots |
| `{prefix}_barplot_superpops.png/.svg` | Standalone continental ancestry barplot |
| `{prefix}_barplot_subpops.png/.svg` | Standalone sub-population barplot (if K=23 available) |

---

## Resource Requirements

| Stage | RAM | Threads | Time (per sample) |
|-------|-----|---------|-------------------|
| fastp | 8 GB | 4 | 5-15 min |
| FastQC | 4 GB | 2 | 5-10 min |
| bwa-mem2 | 48 GB | 16 | 30-90 min |
| HISAT2 | 32 GB | 16 | 20-60 min |
| samtools sort | 32 GB | 8 | 15-30 min |
| MarkDuplicates | 32 GB | 1 | 30-60 min |
| BQSR | 16 GB | 1 | 60-120 min |
| HaplotypeCaller | 16 GB | 4 | 2-12 hours |
| featureCounts | 8 GB | 4 | 5-15 min |
| Ancestry | 80 GB | 4 | 15-60 min |

---

## Reference Setup

All resources live under one directory. Run `scripts/setup_references.sh` to download and index:

```
ngs_resources/
├── genome/                   GRCh38 FASTA + all indices (~35 GB)
├── hisat2_index/             Splice-aware HISAT2 index (~8 GB)
├── annotation/               GENCODE v48 GTF (~1.5 GB)
├── known_sites/              BQSR VCFs from GATK resource bundle (~3.3 GB)
├── intervals/                WES capture BED files (user-provided)
└── ancestry/                 Singularity container (~8 GB)
```

---

## Tools and Versions Used

The exact tool versions this pipeline is configured to load (from the `tools:` block of `config/config.yaml`). This is the reproducibility record — the [Dependencies](#dependencies) table below lists *minimum compatible* versions instead.

| Tool | Version (as configured) | Activation | Stage |
|------|------------------------|-----------|-------|
| Snakemake | 9.17.2 *(run environment; `Snakefile` requires ≥ 8.0)* | conda | Orchestration |
| fastp | 1.0.1 | module | Read QC |
| FastQC | 0.12.1 | module | Read QC |
| MultiQC | 1.20 | module | Read QC |
| bwa-mem2 | 2.3 (conda env `bwa-mem2`, bioconda) | conda | DNA alignment |
| HISAT2 | 2.2.1 | module | RNA alignment |
| SAMtools | 1.21 | module | BAM processing |
| GATK | 4.6.2.0 | module | MarkDuplicates, BQSR, variant calling |
| Subread (featureCounts) | 2.1.1 (conda env `subread`, bioconda) | conda | Quantification |
| Singularity | system binary at `/usr/bin/singularity` *(version not pinned)* | PATH | Ancestry container |
| Python | 3.11 (`python311`) | module | PDF reports (fpdf2, pandas, Pillow) |
| R + Bioconductor | 4.4.2 | module | Ancestry plots (ggplot2, gridExtra; optional readxl) |
| Anaconda | 2023.09 | module | conda environment provisioning |

**Ancestry container (internal tools):** GATK 3.8 + PLINK 1.9 + ADMIXTURE 1.3.0 (as documented in [Stage 5](#stage-5-ancestry-estimation)).

**Reference data:**

| Resource | Version / Source |
|----------|------------------|
| Genome | GRCh38 primary assembly |
| Gene annotation | GENCODE v48 |
| BQSR known sites | GATK resource bundle — dbSNP138, Mills & 1000G gold-standard indels, 1000G Phase 1 high-confidence SNPs, known indels |
| Ancestry reference panel | 1000 Genomes (Phase 3), ancestry-informative markers |

> **Reproducibility note:** `bwa-mem2` and `subread` are activated by conda-env *name*, so the versions above reflect the current envs (`bwa-mem2` 2.3, `subread` 2.1.1 from bioconda) rather than a hard pin in config. Pin them in the env spec if you need them frozen across rebuilds.

---

## Dependencies

| Tool | Version | Stage |
|------|---------|-------|
| Snakemake | >= 8.0 | Orchestration |
| snakemake-executor-plugin-slurm | latest | SLURM job submission (`pip install snakemake-executor-plugin-slurm`) |
| fastp | 0.23+ | Read QC |
| FastQC | 0.12+ | Read QC |
| MultiQC | 1.14+ | Read QC |
| bwa-mem2 | 2.2+ | DNA alignment |
| HISAT2 | 2.2+ | RNA alignment |
| SAMtools | 1.17+ | BAM processing |
| GATK | 4.x | MarkDup, BQSR, variant calling |
| Subread | 2.0+ | featureCounts (RNA) |
| Singularity | 3.x | Ancestry container |
| Python 3 | 3.9+ | Reports (fpdf2, pandas, Pillow) |
| R | 4.x | Ancestry plots (ggplot2, gridExtra; optional: readxl for Excel metadata) |

---

## Authors

Samuel Mwamburi - Yates Lab

---

## Citation

If you find this pipeline useful in your research, please cite it:

> Mwamburi, S. (2026). *NGS Pipeline: a unified Snakemake workflow for RNA-seq, WES, and WGS* (Version 1.0.0) [Computer software]. https://github.com/Kishaz/ngs-pipeline

A machine-readable [`CITATION.cff`](CITATION.cff) is included, so GitHub shows a **"Cite this repository"** button (with APA and BibTeX export) on the repo page.

> There is no associated publication or DOI yet. If you publish work using this pipeline, consider minting a DOI (e.g. via [Zenodo](https://zenodo.org)) and updating the citation above and `CITATION.cff`.

---

## License

Released under the [MIT License](LICENSE) — © 2026 Samuel Mwamburi.
