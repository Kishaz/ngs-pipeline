#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# HISAT2 Alignment Pipeline (RNA-seq)
# Aligns PE or SE reads → sorted, markdup BAM + index + flagstat
#
# Usage:
#   PE:  bash run_hisat2.sh --r1 R1.fq.gz --r2 R2.fq.gz --index grch38 \
#            --sample SampleName --outdir results/alignment \
#            [--threads 16] [--sort-threads 8] [--sort-mem 4G] [--dta]
#
#   SE:  bash run_hisat2.sh --r1 R1.fq.gz --index grch38 \
#            --sample SampleName --outdir results/alignment
###############################################################################

# Defaults
THREADS=16
SORT_THREADS=8
SORT_MEM="4G"
R1=""
R2=""
INDEX=""
SAMPLE=""
OUTDIR=""
DTA="--dta"     # Default: enable for downstream transcript assembly
EXTRA=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --r1 <path>            Trimmed R1 FASTQ (gzipped)
  --index <prefix>       HISAT2 index prefix (e.g., grch38_v48/hisat2_index/grch38)
  --sample <name>        Sample name (used for read group + output naming)
  --outdir <path>        Base output directory

Optional:
  --r2 <path>            Trimmed R2 FASTQ — enables paired-end mode
  --threads <int>        Alignment threads (default: $THREADS)
  --sort-threads <int>   Sorting threads (default: $SORT_THREADS)
  --sort-mem <str>       Memory per sort thread (default: $SORT_MEM)
  --no-dta               Disable --dta flag (for non-transcript-assembly workflows)
  --extra <args>         Additional HISAT2 arguments (quoted string)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)            shift; R1="$1" ;;
        --r2)            shift; R2="$1" ;;
        --index)         shift; INDEX="$1" ;;
        --sample)        shift; SAMPLE="$1" ;;
        --outdir)        shift; OUTDIR="$1" ;;
        --threads)       shift; THREADS="$1" ;;
        --sort-threads)  shift; SORT_THREADS="$1" ;;
        --sort-mem)      shift; SORT_MEM="$1" ;;
        --no-dta)        DTA="" ;;
        --extra)         shift; EXTRA="$1" ;;
        -h|--help)       usage ;;
        *)               echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ -n "$R1" ]]     || { echo "ERROR: --r1 is required"; usage; }
[[ -n "$INDEX" ]]  || { echo "ERROR: --index is required"; usage; }
[[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required"; usage; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; }
[[ -f "$R1" ]]     || { echo "ERROR: R1 file not found: $R1"; exit 1; }

command -v hisat2   >/dev/null || { echo "ERROR: hisat2 not found"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v gatk     >/dev/null || { echo "ERROR: gatk not found"; exit 1; }

# Check index exists (at least one .ht2 file)
if ! ls "${INDEX}"*.ht2 >/dev/null 2>&1 && ! ls "${INDEX}"*.ht2l >/dev/null 2>&1; then
    echo "ERROR: HISAT2 index not found at prefix: $INDEX"
    echo "  Expected files like: ${INDEX}.1.ht2"
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
ALIGN_LOG="${LOG_DIR}/${SAMPLE}.hisat2.log"
ALIGN_SUMMARY="${LOG_DIR}/${SAMPLE}.hisat2_summary.txt"

# Checkpoints
CHECK_ALIGN="${BAM_DIR}/${SAMPLE}.align.done"
CHECK_MARKDUP="${BAM_DIR}/${SAMPLE}.markdup.done"
CHECK_INDEX="${BAM_DIR}/${SAMPLE}.index.done"

echo "========================================="
echo "HISAT2 Alignment Pipeline (RNA-seq)"
echo "========================================="
echo "  Sample:     $SAMPLE"
echo "  Index:      $(basename $INDEX)"
if [[ -n "$R2" ]]; then
    echo "  Mode:       Paired-End"
else
    echo "  Mode:       Single-End"
fi
echo "  Threads:    $THREADS (align) / $SORT_THREADS (sort)"
[[ -n "$DTA" ]] && echo "  DTA:        enabled"
echo "========================================="

# -----------------------------------------------
# STEP 1: ALIGN + SORT
# -----------------------------------------------
if [[ -f "$CHECK_ALIGN" ]]; then
    echo "Alignment already completed — skipping"
else
    echo "Aligning with HISAT2 and sorting..."

    # Build input flags
    if [[ -n "$R2" ]]; then
        [[ -f "$R2" ]] || { echo "ERROR: R2 file not found: $R2"; exit 1; }
        INPUT_FLAGS=(-1 "$R1" -2 "$R2")
    else
        INPUT_FLAGS=(-U "$R1")
    fi

    hisat2 \
        -x "$INDEX" \
        "${INPUT_FLAGS[@]}" \
        $DTA \
        $EXTRA \
        --rg-id "$SAMPLE" \
        --rg "SM:${SAMPLE}" \
        --rg "LB:RNASEQ" \
        --rg "PL:ILLUMINA" \
        --rg "PU:${SAMPLE}" \
        -p "$THREADS" \
        --new-summary \
        --summary-file "$ALIGN_SUMMARY" \
        2> "$ALIGN_LOG" \
    | samtools sort \
        -@ "$SORT_THREADS" \
        -m "$SORT_MEM" \
        -o "$SORTED_BAM" -

    samtools quickcheck "$SORTED_BAM" \
        || { echo "ERROR: BAM integrity check failed for $SAMPLE"; exit 1; }

    touch "$CHECK_ALIGN"
    echo "Alignment complete: $SORTED_BAM"

    # Print alignment summary
    if [[ -f "$ALIGN_SUMMARY" ]]; then
        echo ""
        echo "--- HISAT2 Alignment Summary ---"
        cat "$ALIGN_SUMMARY"
        echo "--------------------------------"
    fi
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
# CLEANUP
# -----------------------------------------------
if [[ -f "$MARKDUP_BAM" ]] && [[ -f "$CHECK_MARKDUP" ]] && [[ -f "$SORTED_BAM" ]]; then
    rm -f "$SORTED_BAM"
    echo "Cleaned up intermediate sorted BAM"
fi

echo ""
echo "========================================="
echo "SAMPLE $SAMPLE COMPLETE"
echo "========================================="
echo "  BAM:            $MARKDUP_BAM"
echo "  Index:          ${MARKDUP_BAM}.bai"
echo "  Flagstat:       $FLAGSTAT"
echo "  Dup Metrics:    $DUP_METRICS"
echo "  HISAT2 Summary: $ALIGN_SUMMARY"
echo "========================================="
