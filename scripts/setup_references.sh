#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# Reference Setup Script for NGS Pipeline
#
# Downloads and indexes all reference resources needed by the pipeline.
# Run this ONCE on your HPC before the first pipeline execution.
#
# Resources are organized under a single RESOURCES_DIR:
#
#   resources/
#   ├── genome/
#   │   ├── GRCh38.primary_assembly.genome.fa        ← Reference FASTA
#   │   ├── GRCh38.primary_assembly.genome.fa.fai    ← samtools index
#   │   ├── GRCh38.primary_assembly.genome.dict      ← GATK dictionary
#   │   ├── GRCh38.primary_assembly.genome.fa.amb     ← bwa-mem2 indices
#   │   ├── GRCh38.primary_assembly.genome.fa.ann
#   │   ├── GRCh38.primary_assembly.genome.fa.pac
#   │   └── GRCh38.primary_assembly.genome.fa.bwt.2bit.64
#   ├── hisat2_index/
#   │   └── grch38.*.ht2                             ← HISAT2 genome index
#   ├── annotation/
#   │   └── gencode.v48.primary_assembly.annotation.gtf
#   ├── known_sites/                                  ← For BQSR (optional)
#   │   ├── dbsnp_146.hg38.vcf.gz
#   │   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#   │   └── 1000G_phase1.snps.high_confidence.hg38.vcf.gz
#   ├── intervals/                                    ← WES capture targets
#   │   └── <your_capture_kit>.bed
#   └── ancestry/
#       └── ancestry-pipeline.sif                     ← Singularity container
#
# Usage:
#   bash setup_references.sh --outdir /projects/yates_lab_hpc/sam/smwambu1/ngs_resources
#
#   # Skip downloads if FASTA already exists:
#   bash setup_references.sh --outdir /path/to/resources --skip-download
#
#   # Only build specific indices:
#   bash setup_references.sh --outdir /path/to/resources --only bwa-mem2
#   bash setup_references.sh --outdir /path/to/resources --only hisat2
#   bash setup_references.sh --outdir /path/to/resources --only all
###############################################################################

OUTDIR=""
SKIP_DOWNLOAD=false
ONLY="all"
THREADS=16
GENCODE_VERSION="48"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --outdir <path>        Base directory for all reference resources

Optional:
  --skip-download        Skip downloading; only build indices from existing files
  --only <target>        Only build specific index: bwa-mem2, hisat2, all (default: all)
  --threads <int>        Threads for indexing (default: $THREADS)
  --gencode-version <N>  GENCODE version (default: $GENCODE_VERSION)

Resource Requirements:
  Disk:   ~60 GB (genome ~3GB + indices ~30GB + annotation ~2GB + known sites ~5GB)
  RAM:    ~80 GB (bwa-mem2 indexing)
  Time:   ~2-4 hours depending on cluster speed
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --outdir)           shift; OUTDIR="$1" ;;
        --skip-download)    SKIP_DOWNLOAD=true ;;
        --only)             shift; ONLY="$1" ;;
        --threads)          shift; THREADS="$1" ;;
        --gencode-version)  shift; GENCODE_VERSION="$1" ;;
        -h|--help)          usage ;;
        *)                  echo "Unknown option: $1"; usage ;;
    esac
    shift
done

[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; }

# Directory structure
GENOME_DIR="${OUTDIR}/genome"
HISAT2_DIR="${OUTDIR}/hisat2_index"
ANNOT_DIR="${OUTDIR}/annotation"
KNOWN_DIR="${OUTDIR}/known_sites"
INTERVAL_DIR="${OUTDIR}/intervals"
ANCESTRY_DIR="${OUTDIR}/ancestry"

mkdir -p "$GENOME_DIR" "$HISAT2_DIR" "$ANNOT_DIR" "$KNOWN_DIR" "$INTERVAL_DIR" "$ANCESTRY_DIR"

REF_FA="${GENOME_DIR}/GRCh38.primary_assembly.genome.fa"
GTF="${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf"

echo "================================================================"
echo "  NGS Pipeline — Reference Setup"
echo "================================================================"
echo "  Output:    $OUTDIR"
echo "  Threads:   $THREADS"
echo "  GENCODE:   v${GENCODE_VERSION}"
echo "  Download:  $([ "$SKIP_DOWNLOAD" = true ] && echo 'skipped' || echo 'enabled')"
echo "  Target:    $ONLY"
echo "================================================================"
echo ""

###############################################################################
# STEP 1: Download reference genome (GRCh38 primary assembly)
###############################################################################
if [[ "$SKIP_DOWNLOAD" = false ]]; then

    echo "--- Step 1: Reference Genome ---"
    if [[ -f "$REF_FA" ]]; then
        echo "  FASTA already exists — skipping download"
    else
        echo "  Downloading GRCh38 primary assembly from GENCODE..."
        wget -q --show-progress -O "${REF_FA}.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa.gz"

        echo "  Decompressing..."
        gunzip "${REF_FA}.gz"
        echo "  Done: $REF_FA"
    fi
    echo ""

    ###########################################################################
    # STEP 2: Download gene annotation (GTF)
    ###########################################################################
    echo "--- Step 2: Gene Annotation (GTF) ---"
    if [[ -f "$GTF" ]]; then
        echo "  GTF already exists — skipping download"
    else
        echo "  Downloading GENCODE v${GENCODE_VERSION} annotation..."
        wget -q --show-progress -O "${GTF}.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf.gz"

        echo "  Decompressing..."
        gunzip "${GTF}.gz"
        echo "  Done: $GTF"
    fi
    echo ""

    ###########################################################################
    # STEP 3: Download known sites for BQSR (optional but recommended)
    ###########################################################################
    echo "--- Step 3: Known Sites (BQSR) ---"
    echo "  These files are from the GATK resource bundle."
    echo "  If you have them already, copy them to: $KNOWN_DIR"
    echo ""
    echo "  Expected files (download manually from GATK bundle if needed):"
    echo "    - dbsnp_146.hg38.vcf.gz (+.tbi)"
    echo "    - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz (+.tbi)"
    echo "    - 1000G_phase1.snps.high_confidence.hg38.vcf.gz (+.tbi)"
    echo ""
    echo "  GATK resource bundle:"
    echo "    gs://genomics-public-data/resources/broad/hg38/v0/"
    echo "    https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0"
    echo ""

fi

###############################################################################
# STEP 4: Build indices
###############################################################################
[[ -f "$REF_FA" ]] || { echo "ERROR: Reference FASTA not found: $REF_FA"; exit 1; }

# --- samtools faidx (always needed) ---
echo "--- Step 4a: samtools faidx ---"
if [[ -f "${REF_FA}.fai" ]]; then
    echo "  .fai index already exists — skipping"
else
    echo "  Building FASTA index..."
    samtools faidx "$REF_FA"
    echo "  Done: ${REF_FA}.fai"
fi
echo ""

# --- GATK sequence dictionary (always needed) ---
echo "--- Step 4b: GATK sequence dictionary ---"
DICT="${GENOME_DIR}/GRCh38.primary_assembly.genome.dict"
if [[ -f "$DICT" ]]; then
    echo "  .dict already exists — skipping"
else
    echo "  Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R "$REF_FA"
    echo "  Done: $DICT"
fi
echo ""

# --- bwa-mem2 index (for DNA alignment) ---
if [[ "$ONLY" == "all" || "$ONLY" == "bwa-mem2" ]]; then
    echo "--- Step 4c: bwa-mem2 index ---"
    echo "  WARNING: Requires ~80 GB RAM and ~1-2 hours"
    if [[ -f "${REF_FA}.bwt.2bit.64" ]] || [[ -f "${REF_FA}.bwt.2bit" ]]; then
        echo "  bwa-mem2 index already exists — skipping"
    else
        echo "  Building bwa-mem2 index..."
        bwa-mem2 index "$REF_FA"
        echo "  Done"
    fi
    echo ""
fi

# --- HISAT2 index (for RNA alignment) ---
if [[ "$ONLY" == "all" || "$ONLY" == "hisat2" ]]; then
    echo "--- Step 4d: HISAT2 index ---"
    echo "  WARNING: Requires ~200 GB RAM and ~2-3 hours"
    HISAT2_PREFIX="${HISAT2_DIR}/grch38"
    if ls "${HISAT2_PREFIX}"*.ht2 >/dev/null 2>&1 || ls "${HISAT2_PREFIX}"*.ht2l >/dev/null 2>&1; then
        echo "  HISAT2 index already exists — skipping"
    else
        echo "  Building HISAT2 index..."
        echo "  (For faster builds, use hisat2_extract_splice_sites.py and"
        echo "   hisat2_extract_exons.py with the GTF for splice-aware indexing)"
        echo ""

        if [[ -f "$GTF" ]]; then
            echo "  Extracting splice sites and exons from GTF..."
            hisat2_extract_splice_sites.py "$GTF" > "${HISAT2_DIR}/splice_sites.tsv"
            hisat2_extract_exons.py "$GTF" > "${HISAT2_DIR}/exons.tsv"

            echo "  Building splice-aware index..."
            hisat2-build \
                -p "$THREADS" \
                --ss "${HISAT2_DIR}/splice_sites.tsv" \
                --exon "${HISAT2_DIR}/exons.tsv" \
                "$REF_FA" \
                "$HISAT2_PREFIX"
        else
            echo "  GTF not found — building basic index (no splice sites)..."
            hisat2-build -p "$THREADS" "$REF_FA" "$HISAT2_PREFIX"
        fi
        echo "  Done: ${HISAT2_PREFIX}.*.ht2"
    fi
    echo ""
fi

###############################################################################
# STEP 4b: RNA-seq QC references (BED12 / refFlat / rRNA intervals)
# Built from the GTF with a self-contained Python converter (no UCSC tools,
# no network). Used by RSeQC (strandedness, read distribution) and Picard
# CollectRnaSeqMetrics. Config keys: rnaseq_qc.{bed12,refflat,rrna_intervals}.
###############################################################################
echo "--- Step 4b: RNA-seq QC references ---"
BED12="${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.genes.bed12"
REFFLAT="${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.refFlat.txt"
RRNA="${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.rRNA.interval_list"
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$BED12" && -f "$REFFLAT" && -f "$RRNA" ]]; then
    echo "  RNA-seq QC references already exist — skipping"
elif [[ -f "$GTF" ]]; then
    python3 "${_SCRIPT_DIR}/gtf_to_rnaseq_refs.py" \
        --gtf "$GTF" \
        --dict "$DICT" \
        --bed12 "$BED12" \
        --refflat "$REFFLAT" \
        --rrna "$RRNA"
    echo "  Done: $BED12 / $REFFLAT / $RRNA"
else
    echo "  GTF not found — skipping RNA-seq QC references"
fi
echo ""

###############################################################################
# STEP 5: Verify everything
###############################################################################
echo "================================================================"
echo "  Reference Setup — Verification"
echo "================================================================"

check_file() {
    local label="$1"
    local path="$2"
    if [[ -f "$path" ]]; then
        local size=$(du -sh "$path" 2>/dev/null | cut -f1)
        printf "  %-45s %8s   OK\n" "$label" "$size"
    elif [[ -d "$path" ]]; then
        printf "  %-45s            OK (dir)\n" "$label"
    else
        printf "  %-45s            MISSING\n" "$label"
    fi
}

echo ""
echo "  Genome:"
check_file "Reference FASTA"        "$REF_FA"
check_file "FASTA index (.fai)"     "${REF_FA}.fai"
check_file "Sequence dict (.dict)"  "$DICT"
echo ""

echo "  bwa-mem2 indices:"
check_file ".amb"                   "${REF_FA}.amb"
check_file ".ann"                   "${REF_FA}.ann"
check_file ".pac"                   "${REF_FA}.pac"
check_file ".bwt.2bit.64"          "${REF_FA}.bwt.2bit.64"
echo ""

echo "  HISAT2 index:"
HISAT2_PREFIX="${HISAT2_DIR}/grch38"
if ls "${HISAT2_PREFIX}"*.ht2 >/dev/null 2>&1; then
    N_HT2=$(ls "${HISAT2_PREFIX}"*.ht2 2>/dev/null | wc -l | tr -d ' ')
    printf "  %-45s            OK (%s files)\n" "HISAT2 index files" "$N_HT2"
elif ls "${HISAT2_PREFIX}"*.ht2l >/dev/null 2>&1; then
    N_HT2=$(ls "${HISAT2_PREFIX}"*.ht2l 2>/dev/null | wc -l | tr -d ' ')
    printf "  %-45s            OK (%s files, large index)\n" "HISAT2 index files" "$N_HT2"
else
    printf "  %-45s            MISSING\n" "HISAT2 index files"
fi
echo ""

echo "  Annotation:"
check_file "GENCODE GTF"            "$GTF"
check_file "RNA QC BED12"          "${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.genes.bed12"
check_file "RNA QC refFlat"        "${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.refFlat.txt"
check_file "RNA QC rRNA intervals" "${ANNOT_DIR}/gencode.v${GENCODE_VERSION}.rRNA.interval_list"
echo ""

echo "  Known sites (BQSR):"
check_file "dbSNP"                  "${KNOWN_DIR}/dbsnp_146.hg38.vcf.gz"
check_file "Mills & 1000G indels"   "${KNOWN_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
check_file "1000G high-conf SNPs"   "${KNOWN_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
echo ""

echo "  WES intervals:"
BED_FILES=( "${INTERVAL_DIR}"/*.bed )
if [[ ${#BED_FILES[@]} -gt 0 ]] && [[ -f "${BED_FILES[0]}" ]]; then
    for bed in "${BED_FILES[@]}"; do
        check_file "$(basename $bed)" "$bed"
    done
else
    printf "  %-45s            NONE (add your capture kit BED)\n" "Target BED files"
fi
echo ""

echo "  Ancestry container:"
SIF_FILES=( "${ANCESTRY_DIR}"/*.sif )
if [[ ${#SIF_FILES[@]} -gt 0 ]] && [[ -f "${SIF_FILES[0]}" ]]; then
    for sif in "${SIF_FILES[@]}"; do
        check_file "$(basename $sif)" "$sif"
    done
else
    printf "  %-45s            NONE\n" "Singularity .sif image"
fi

echo ""
echo "================================================================"
echo ""
echo "  Update your pipeline config with these paths:"
echo ""
echo "    ref:"
echo "      fasta: \"${REF_FA}\""
echo ""
echo "    hisat2_index: \"${HISAT2_DIR}/grch38\""
echo "    bwa_mem2_index: \"${REF_FA}\""
echo "    gtf: \"${GTF}\""
echo ""
echo "    ancestry:"
echo "      container: \"${ANCESTRY_DIR}/ancestry-pipeline.sif\""
echo ""
echo "    variant_calling:"
echo "      haplotypecaller:"
echo "        intervals: \"${INTERVAL_DIR}/<your_capture_kit>.bed\"  # WES only"
echo "      genomicsdb:"
echo "        intervals: \"${INTERVAL_DIR}/<your_intervals>.list\"   # Required for joint calling"
echo ""
echo "================================================================"
