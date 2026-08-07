#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# fastp — Adapter removal + quality trimming
# Supports: Paired-End and Single-End reads
#
# Usage:
#   PE:  bash run_fastp.sh --r1 R1.fq.gz --r2 R2.fq.gz --out-r1 trimmed_R1.fq.gz \
#            --out-r2 trimmed_R2.fq.gz --json report.json --html report.html \
#            [--qual 20] [--minlen 50] [--threads 4]
#
#   SE:  bash run_fastp.sh --r1 R1.fq.gz --out-r1 trimmed_R1.fq.gz \
#            --json report.json --html report.html \
#            [--qual 20] [--minlen 50] [--threads 4]
###############################################################################

# Defaults
QUAL=20
MINLEN=50
THREADS=4
R1=""
R2=""
OUT_R1=""
OUT_R2=""
JSON=""
HTML=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --r1 <path>        Input R1 FASTQ (gzipped)
  --out-r1 <path>    Output trimmed R1 FASTQ
  --json <path>      Output fastp JSON report
  --html <path>      Output fastp HTML report

Optional (PE mode):
  --r2 <path>        Input R2 FASTQ (gzipped) — enables paired-end mode
  --out-r2 <path>    Output trimmed R2 FASTQ (required if --r2 given)

Parameters:
  --qual <int>       Minimum base quality (default: $QUAL)
  --minlen <int>     Minimum read length after trimming (default: $MINLEN)
  --threads <int>    Number of threads (default: $THREADS)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)        shift; R1="$1" ;;
        --r2)        shift; R2="$1" ;;
        --out-r1)    shift; OUT_R1="$1" ;;
        --out-r2)    shift; OUT_R2="$1" ;;
        --json)      shift; JSON="$1" ;;
        --html)      shift; HTML="$1" ;;
        --qual)      shift; QUAL="$1" ;;
        --minlen)    shift; MINLEN="$1" ;;
        --threads)   shift; THREADS="$1" ;;
        -h|--help)   usage ;;
        *)           echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ -n "$R1" ]]     || { echo "ERROR: --r1 is required"; usage; }
[[ -n "$OUT_R1" ]] || { echo "ERROR: --out-r1 is required"; usage; }
[[ -n "$JSON" ]]   || { echo "ERROR: --json is required"; usage; }
[[ -n "$HTML" ]]   || { echo "ERROR: --html is required"; usage; }
[[ -f "$R1" ]]     || { echo "ERROR: R1 file not found: $R1"; exit 1; }

command -v fastp >/dev/null || { echo "ERROR: fastp not found in PATH"; exit 1; }

# Create output directories
mkdir -p "$(dirname "$OUT_R1")" "$(dirname "$JSON")"

# Determine mode
if [[ -n "$R2" ]]; then
    # ---- Paired-End ----
    [[ -f "$R2" ]]     || { echo "ERROR: R2 file not found: $R2"; exit 1; }
    [[ -n "$OUT_R2" ]] || { echo "ERROR: --out-r2 required for paired-end mode"; usage; }

    SAMPLE=$(basename "$R1" | sed 's/_R1.*$//')
    echo "-----------------------------------------"
    echo "fastp PE | Sample: ${SAMPLE}"
    echo "  R1: $R1"
    echo "  R2: $R2"
    echo "  Quality >= $QUAL | Min length >= $MINLEN | Threads: $THREADS"
    echo "-----------------------------------------"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT_R1" \
        -O "$OUT_R2" \
        --json "$JSON" \
        --html "$HTML" \
        --detect_adapter_for_pe \
        --qualified_quality_phred "$QUAL" \
        --length_required "$MINLEN" \
        --thread "$THREADS"

    echo "Trimmed PE reads:"
    echo "  R1: $OUT_R1"
    echo "  R2: $OUT_R2"
else
    # ---- Single-End ----
    SAMPLE=$(basename "$R1" | sed 's/_R1.*$//')
    echo "-----------------------------------------"
    echo "fastp SE | Sample: ${SAMPLE}"
    echo "  R1: $R1"
    echo "  Quality >= $QUAL | Min length >= $MINLEN | Threads: $THREADS"
    echo "-----------------------------------------"

    fastp \
        -i "$R1" \
        -o "$OUT_R1" \
        --json "$JSON" \
        --html "$HTML" \
        --qualified_quality_phred "$QUAL" \
        --length_required "$MINLEN" \
        --thread "$THREADS"

    echo "Trimmed SE reads:"
    echo "  R1: $OUT_R1"
fi

echo "Reports: $JSON, $HTML"
echo "fastp complete."
