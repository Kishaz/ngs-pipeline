#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# FastQC — Per-file quality assessment
#
# Usage:
#   bash run_fastqc.sh --input R1.fastq.gz --outdir qc/fastqc_raw [--threads 2]
#
# Accepts one or more --input flags for batch processing.
###############################################################################

THREADS=2
OUTDIR=""
INPUTS=()

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --input <path>     Input FASTQ file(s) — can be specified multiple times
  --outdir <path>    Output directory for FastQC reports

Optional:
  --threads <int>    Number of threads (default: $THREADS)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --input)     shift; INPUTS+=("$1") ;;
        --outdir)    shift; OUTDIR="$1" ;;
        --threads)   shift; THREADS="$1" ;;
        -h|--help)   usage ;;
        *)           echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "ERROR: at least one --input is required"; usage; }
[[ -n "$OUTDIR" ]]        || { echo "ERROR: --outdir is required"; usage; }

command -v fastqc >/dev/null || { echo "ERROR: fastqc not found in PATH"; exit 1; }

mkdir -p "$OUTDIR"

# Verify all inputs exist
for f in "${INPUTS[@]}"; do
    [[ -f "$f" ]] || { echo "ERROR: Input file not found: $f"; exit 1; }
done

echo "-----------------------------------------"
echo "FastQC | ${#INPUTS[@]} file(s) | Threads: $THREADS"
echo "  Output: $OUTDIR"
echo "-----------------------------------------"

fastqc \
    -t "$THREADS" \
    -o "$OUTDIR" \
    --noextract \
    "${INPUTS[@]}"

echo ""
echo "FastQC reports written to: $OUTDIR"
for f in "${INPUTS[@]}"; do
    base=$(basename "$f" .fastq.gz)
    echo "  ${base}_fastqc.html"
done
echo "FastQC complete."
