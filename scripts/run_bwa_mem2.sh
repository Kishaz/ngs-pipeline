#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# bwa-mem2 Alignment Pipeline (DNA: WES / WGS)
# Aligns PE or SE reads → sorted, markdup BAM + index + flagstat
#
# Usage:
#   PE:  bash run_bwa_mem2.sh --r1 R1.fq.gz --r2 R2.fq.gz --ref ref.fa \
#            --sample SampleName --outdir results/alignment \
#            --seq-type WES [--threads 16] [--sort-threads 8] [--sort-mem 4G]
#
#   SE:  bash run_bwa_mem2.sh --r1 R1.fq.gz --ref ref.fa \
#            --sample SampleName --outdir results/alignment \
#            --seq-type WGS
###############################################################################

# Defaults
THREADS=16
SORT_THREADS=8
SORT_MEM="4G"
R1=""
R2=""
REF=""
SAMPLE=""
OUTDIR=""
SEQ_TYPE="WES"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --r1 <path>            Trimmed R1 FASTQ (gzipped)
  --ref <path>           Reference FASTA (bwa-mem2 indexed)
  --sample <name>        Sample name (used for read group + output naming)
  --outdir <path>        Base output directory

Optional:
  --r2 <path>            Trimmed R2 FASTQ — enables paired-end mode
  --seq-type <WES|WGS>   Sequencing type for read group LB tag (default: $SEQ_TYPE)
  --threads <int>        Alignment threads (default: $THREADS)
  --sort-threads <int>   Sorting threads (default: $SORT_THREADS)
  --sort-mem <str>       Memory per sort thread (default: $SORT_MEM)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)            shift; R1="$1" ;;
        --r2)            shift; R2="$1" ;;
        --ref)           shift; REF="$1" ;;
        --sample)        shift; SAMPLE="$1" ;;
        --outdir)        shift; OUTDIR="$1" ;;
        --seq-type)      shift; SEQ_TYPE="$1" ;;
        --threads)       shift; THREADS="$1" ;;
        --sort-threads)  shift; SORT_THREADS="$1" ;;
        --sort-mem)      shift; SORT_MEM="$1" ;;
        -h|--help)       usage ;;
        *)               echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ -n "$R1" ]]     || { echo "ERROR: --r1 is required"; usage; }
[[ -n "$REF" ]]    || { echo "ERROR: --ref is required"; usage; }
[[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required"; usage; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; }
[[ -f "$R1" ]]     || { echo "ERROR: R1 file not found: $R1"; exit 1; }
[[ -f "$REF" ]]    || { echo "ERROR: Reference not found: $REF"; exit 1; }

command -v bwa-mem2 >/dev/null || { echo "ERROR: bwa-mem2 not found"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v gatk     >/dev/null || { echo "ERROR: gatk not found"; exit 1; }

# Reference index checks
[[ -f "${REF}.fai" ]] || { echo "ERROR: Missing .fai index for reference"; exit 1; }
for ext in amb ann pac; do
    [[ -f "${REF}.${ext}" ]] || { echo "ERROR: Missing bwa index (.${ext})"; exit 1; }
done
if [[ ! -f "${REF}.bwt.2bit.64" ]] && [[ ! -f "${REF}.bwt.2bit" ]]; then
    echo "ERROR: Missing bwa-mem2 BWT index (.bwt.2bit or .bwt.2bit.64)"
    exit 1
fi

# Directories
BAM_DIR="${OUTDIR}/bam"
METRICS_DIR="${OUTDIR}/metrics"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "$BAM_DIR" "$METRICS_DIR" "$LOG_DIR"

# Output paths
SORTED_BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
MARKDUP_BAM="${BAM_DIR}/${SAMPLE}.sorted.markdup.bam"
DUP_METRICS="${METRICS_DIR}/${SAMPLE}.dup_metrics.txt"
FLAGSTAT="${METRICS_DIR}/${SAMPLE}.flagstat.txt"
ALIGN_LOG="${LOG_DIR}/${SAMPLE}.bwa_mem2.log"

# Checkpoints
CHECK_ALIGN="${BAM_DIR}/${SAMPLE}.align.done"
CHECK_MARKDUP="${BAM_DIR}/${SAMPLE}.markdup.done"
CHECK_INDEX="${BAM_DIR}/${SAMPLE}.index.done"

# Read group
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SEQ_TYPE}\tPL:ILLUMINA\tPU:${SAMPLE}"

echo "========================================="
echo "bwa-mem2 Alignment Pipeline"
echo "========================================="
echo "  Sample:     $SAMPLE"
echo "  Seq Type:   $SEQ_TYPE"
echo "  Reference:  $(basename $REF)"
if [[ -n "$R2" ]]; then
    echo "  Mode:       Paired-End"
else
    echo "  Mode:       Single-End"
fi
echo "  Threads:    $THREADS (align) / $SORT_THREADS (sort)"
echo "========================================="

# -----------------------------------------------
# STEP 1: ALIGN + SORT
# -----------------------------------------------
if [[ -f "$CHECK_ALIGN" ]]; then
    echo "Alignment already completed — skipping"
else
    echo "Aligning and sorting..."

    if [[ -n "$R2" ]]; then
        bwa-mem2 mem \
            -Y \
            -t "$THREADS" \
            -R "$RG" \
            "$REF" "$R1" "$R2" 2> "$ALIGN_LOG" \
        | samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEM" \
            -o "$SORTED_BAM" -
    else
        bwa-mem2 mem \
            -Y \
            -t "$THREADS" \
            -R "$RG" \
            "$REF" "$R1" 2> "$ALIGN_LOG" \
        | samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEM" \
            -o "$SORTED_BAM" -
    fi

    samtools quickcheck "$SORTED_BAM" \
        || { echo "ERROR: BAM integrity check failed for $SAMPLE"; exit 1; }

    touch "$CHECK_ALIGN"
    echo "Alignment complete: $SORTED_BAM"
fi

# -----------------------------------------------
# STEP 2: MARK DUPLICATES
# -----------------------------------------------
if [[ -f "$CHECK_MARKDUP" ]]; then
    echo "MarkDuplicates already completed — skipping"
else
    echo "Marking duplicates..."

    gatk MarkDuplicates \
        -I "$SORTED_BAM" \
        -O "$MARKDUP_BAM" \
        -M "$DUP_METRICS" \
        --CREATE_INDEX false \
        --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY LENIENT

    touch "$CHECK_MARKDUP"
    echo "MarkDuplicates complete: $MARKDUP_BAM"
fi

# -----------------------------------------------
# STEP 3: INDEX
# -----------------------------------------------
if [[ -f "$CHECK_INDEX" ]]; then
    echo "BAM index already exists — skipping"
else
    echo "Indexing BAM..."
    samtools index -@ "$SORT_THREADS" "$MARKDUP_BAM"
    touch "$CHECK_INDEX"
fi

# -----------------------------------------------
# STEP 4: FLAGSTAT
# -----------------------------------------------
if [[ ! -f "$FLAGSTAT" ]]; then
    echo "Computing flagstat..."
    samtools flagstat "$MARKDUP_BAM" > "$FLAGSTAT"
fi

# -----------------------------------------------
# CLEANUP — remove intermediate sorted BAM
# -----------------------------------------------
if [[ -f "$MARKDUP_BAM" ]] && [[ -f "$CHECK_MARKDUP" ]] && [[ -f "$SORTED_BAM" ]]; then
    rm -f "$SORTED_BAM"
    echo "Cleaned up intermediate sorted BAM"
fi

echo ""
echo "========================================="
echo "SAMPLE $SAMPLE COMPLETE"
echo "========================================="
echo "  BAM:         $MARKDUP_BAM"
echo "  Index:       ${MARKDUP_BAM}.bai"
echo "  Flagstat:    $FLAGSTAT"
echo "  Dup Metrics: $DUP_METRICS"
echo "========================================="
