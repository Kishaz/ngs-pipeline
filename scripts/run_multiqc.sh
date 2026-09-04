#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# MultiQC — Aggregate QC reports across samples
#
# Usage:
#   bash run_multiqc.sh --indir results/qc --outdir results/qc/multiqc \
#       [--title "Project QC"] [--config multiqc_config.yaml]
###############################################################################

INDIRS=()
OUTDIR=""
TITLE=""
CONFIG=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --indir <path>     Directory to scan for QC reports. Repeatable — pass it
                     multiple times to aggregate several trees (e.g. qc,
                     alignment/metrics, counts) into one report.
  --outdir <path>    Output directory for MultiQC report

Optional:
  --title <string>   Report title
  --config <path>    Custom MultiQC config YAML
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --indir)     shift; INDIRS+=("$1") ;;
        --outdir)    shift; OUTDIR="$1" ;;
        --title)     shift; TITLE="$1" ;;
        --config)    shift; CONFIG="$1" ;;
        -h|--help)   usage ;;
        *)           echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ ${#INDIRS[@]} -gt 0 ]] || { echo "ERROR: at least one --indir is required"; usage; }
[[ -n "$OUTDIR" ]]        || { echo "ERROR: --outdir is required"; usage; }

command -v multiqc >/dev/null || { echo "ERROR: multiqc not found in PATH"; exit 1; }

# Keep only directories that actually exist (a stage may not have produced some).
SCAN=()
for d in "${INDIRS[@]}"; do
    if [[ -d "$d" ]]; then SCAN+=("$d"); else echo "NOTE: skipping missing input dir: $d"; fi
done
[[ ${#SCAN[@]} -gt 0 ]] || { echo "ERROR: none of the input directories exist"; exit 1; }

mkdir -p "$OUTDIR"

# Build command as array. --ignore keeps MultiQC from walking large BAM trees.
CMD=(multiqc "${SCAN[@]}" -o "$OUTDIR" --force --ignore "*/bam/*")
[[ -n "$TITLE" ]]  && CMD+=(--title "$TITLE")
[[ -n "$CONFIG" ]] && CMD+=(-c "$CONFIG")

echo "-----------------------------------------"
echo "MultiQC"
echo "  Inputs: ${SCAN[*]}"
echo "  Output: $OUTDIR"
[[ -n "$TITLE" ]]  && echo "  Title:  $TITLE"
[[ -n "$CONFIG" ]] && echo "  Config: $CONFIG"
echo "-----------------------------------------"

"${CMD[@]}"

echo ""
echo "MultiQC report: ${OUTDIR}/multiqc_report.html"
echo "MultiQC complete."
