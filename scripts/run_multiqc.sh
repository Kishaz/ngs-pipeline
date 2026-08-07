#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# MultiQC — Aggregate QC reports across samples
#
# Usage:
#   bash run_multiqc.sh --indir results/qc --outdir results/qc/multiqc \
#       [--title "Project QC"] [--config multiqc_config.yaml]
###############################################################################

INDIR=""
OUTDIR=""
TITLE=""
CONFIG=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --indir <path>     Directory containing FastQC/fastp reports to aggregate
  --outdir <path>    Output directory for MultiQC report

Optional:
  --title <string>   Report title
  --config <path>    Custom MultiQC config YAML
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --indir)     shift; INDIR="$1" ;;
        --outdir)    shift; OUTDIR="$1" ;;
        --title)     shift; TITLE="$1" ;;
        --config)    shift; CONFIG="$1" ;;
        -h|--help)   usage ;;
        *)           echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ -n "$INDIR" ]]  || { echo "ERROR: --indir is required"; usage; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; }
[[ -d "$INDIR" ]]  || { echo "ERROR: Input directory not found: $INDIR"; exit 1; }

command -v multiqc >/dev/null || { echo "ERROR: multiqc not found in PATH"; exit 1; }

mkdir -p "$OUTDIR"

# Build command as array
CMD=(multiqc "$INDIR" -o "$OUTDIR" --force)
[[ -n "$TITLE" ]]  && CMD+=(--title "$TITLE")
[[ -n "$CONFIG" ]] && CMD+=(-c "$CONFIG")

echo "-----------------------------------------"
echo "MultiQC"
echo "  Input:  $INDIR"
echo "  Output: $OUTDIR"
[[ -n "$TITLE" ]]  && echo "  Title:  $TITLE"
[[ -n "$CONFIG" ]] && echo "  Config: $CONFIG"
echo "-----------------------------------------"

"${CMD[@]}"

echo ""
echo "MultiQC report: ${OUTDIR}/multiqc_report.html"
echo "MultiQC complete."
