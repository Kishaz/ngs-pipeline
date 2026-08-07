#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# GATK Base Quality Score Recalibration (BQSR) — DNA only
#
# Two-step process:
#   1. BaseRecalibrator  → models systematic base quality errors using known sites
#   2. ApplyBQSR         → corrects base qualities in the BAM
#
# Usage:
#   bash run_bqsr.sh \
#       --bam sample.sorted.markdup.bam \
#       --ref GRCh38.fa \
#       --sample SampleName \
#       --outdir results/alignment \
#       --known-sites dbsnp138.vcf.gz \
#       --known-sites Mills_indels.vcf.gz \
#       --known-sites known_indels.vcf.gz \
#       [--intervals targets.bed]
###############################################################################

BAM=""
REF=""
SAMPLE=""
OUTDIR=""
KNOWN_SITES=()
INTERVALS=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --bam <path>             Input BAM (coordinate-sorted, markdup, indexed)
  --ref <path>             Reference FASTA (with .fai and .dict)
  --sample <name>          Sample name
  --outdir <path>          Output directory
  --known-sites <path>     Known variant sites VCF (can be specified multiple times)

Optional:
  --intervals <path>       Restrict to genomic regions (BED file for WES)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)          shift; BAM="$1" ;;
        --ref)          shift; REF="$1" ;;
        --sample)       shift; SAMPLE="$1" ;;
        --outdir)       shift; OUTDIR="$1" ;;
        --known-sites)  shift; KNOWN_SITES+=("$1") ;;
        --intervals)    shift; INTERVALS="$1" ;;
        -h|--help)      usage ;;
        *)              echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ -n "$BAM" ]]                   || { echo "ERROR: --bam is required"; usage; }
[[ -n "$REF" ]]                   || { echo "ERROR: --ref is required"; usage; }
[[ -n "$SAMPLE" ]]                || { echo "ERROR: --sample is required"; usage; }
[[ -n "$OUTDIR" ]]                || { echo "ERROR: --outdir is required"; usage; }
[[ ${#KNOWN_SITES[@]} -gt 0 ]]   || { echo "ERROR: at least one --known-sites is required"; usage; }
[[ -f "$BAM" ]]                   || { echo "ERROR: BAM not found: $BAM"; exit 1; }
[[ -f "$REF" ]]                   || { echo "ERROR: Reference not found: $REF"; exit 1; }

command -v gatk     >/dev/null || { echo "ERROR: gatk not found"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }

# Check BAM index
if [[ ! -f "${BAM}.bai" ]] && [[ ! -f "${BAM%.*}.bai" ]]; then
    echo "ERROR: BAM index not found. Run: samtools index $BAM"
    exit 1
fi

# Check known sites exist
for ks in "${KNOWN_SITES[@]}"; do
    [[ -f "$ks" ]] || { echo "ERROR: Known sites file not found: $ks"; exit 1; }
done

# Directories
BQSR_DIR="${OUTDIR}/bqsr"
BAM_DIR="${OUTDIR}/bam"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "$BQSR_DIR" "$BAM_DIR" "$LOG_DIR"

# Output paths
RECAL_TABLE="${BQSR_DIR}/${SAMPLE}.recal_data.table"
RECAL_BAM="${BAM_DIR}/${SAMPLE}.sorted.markdup.recal.bam"
LOG_BR="${LOG_DIR}/${SAMPLE}.base_recalibrator.log"
LOG_APPLY="${LOG_DIR}/${SAMPLE}.apply_bqsr.log"

# Checkpoints
CHECK_BR="${BQSR_DIR}/${SAMPLE}.br.done"
CHECK_APPLY="${BQSR_DIR}/${SAMPLE}.apply.done"

# Build known-sites args
KS_ARGS=""
for ks in "${KNOWN_SITES[@]}"; do
    KS_ARGS+=" --known-sites $ks"
done

# Build intervals arg
INTERVAL_ARG=""
if [[ -n "$INTERVALS" ]]; then
    [[ -f "$INTERVALS" ]] || { echo "ERROR: Intervals file not found: $INTERVALS"; exit 1; }
    INTERVAL_ARG="-L $INTERVALS"
fi

echo "========================================="
echo "GATK BQSR Pipeline"
echo "========================================="
echo "  Sample:      $SAMPLE"
echo "  Input BAM:   $(basename $BAM)"
echo "  Reference:   $(basename $REF)"
echo "  Known sites: ${#KNOWN_SITES[@]} file(s)"
for ks in "${KNOWN_SITES[@]}"; do
    echo "    - $(basename $ks)"
done
echo "  Intervals:   ${INTERVALS:-whole genome}"
echo "========================================="

# -----------------------------------------------
# STEP 1: BaseRecalibrator
# -----------------------------------------------
if [[ -f "$CHECK_BR" ]]; then
    echo "BaseRecalibrator already completed — skipping"
else
    echo ""
    echo "Step 1/2: BaseRecalibrator..."
    echo "  Building recalibration model from known variant sites"

    gatk BaseRecalibrator \
        -R "$REF" \
        -I "$BAM" \
        $KS_ARGS \
        $INTERVAL_ARG \
        -O "$RECAL_TABLE" \
        2> "$LOG_BR"

    [[ -f "$RECAL_TABLE" ]] || { echo "ERROR: Recalibration table not created"; exit 1; }

    touch "$CHECK_BR"
    echo "  Recalibration table: $RECAL_TABLE"
fi

# -----------------------------------------------
# STEP 2: ApplyBQSR
# -----------------------------------------------
if [[ -f "$CHECK_APPLY" ]]; then
    echo "ApplyBQSR already completed — skipping"
else
    echo ""
    echo "Step 2/2: ApplyBQSR..."
    echo "  Applying base quality recalibration"

    gatk ApplyBQSR \
        -R "$REF" \
        -I "$BAM" \
        --bqsr-recal-file "$RECAL_TABLE" \
        -O "$RECAL_BAM" \
        2> "$LOG_APPLY"

    [[ -f "$RECAL_BAM" ]] || { echo "ERROR: Recalibrated BAM not created"; exit 1; }

    touch "$CHECK_APPLY"
    echo "  Recalibrated BAM: $RECAL_BAM"
fi

# -----------------------------------------------
# STEP 3: Index recalibrated BAM
# -----------------------------------------------
if [[ ! -f "${RECAL_BAM}.bai" ]]; then
    echo ""
    echo "Indexing recalibrated BAM..."
    samtools index "$RECAL_BAM"
fi

echo ""
echo "========================================="
echo "BQSR COMPLETE — $SAMPLE"
echo "========================================="
echo "  Recal table: $RECAL_TABLE"
echo "  Recal BAM:   $RECAL_BAM"
echo "  Recal index: ${RECAL_BAM}.bai"
echo "========================================="
