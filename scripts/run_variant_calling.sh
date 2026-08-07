#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# GATK Variant Calling Pipeline (DNA: WES / WGS)
#
# Per-sample:  HaplotypeCaller → per-sample gVCF
# Joint:       GenomicsDBImport → GenotypeGVCFs → cohort VCF
#
# Modes:
#   --mode single    Run HaplotypeCaller on one sample (produces .g.vcf.gz)
#   --mode joint     Run GenomicsDBImport + GenotypeGVCFs on all gVCFs
#
# Usage:
#   # Per-sample gVCF
#   bash run_variant_calling.sh --mode single \
#       --bam sample.sorted.markdup.bam \
#       --ref GRCh38.fa \
#       --sample SampleName \
#       --outdir results/variants \
#       [--intervals targets.bed] [--threads 4]
#
#   # Joint genotyping
#   bash run_variant_calling.sh --mode joint \
#       --gvcf-dir results/variants/gvcf \
#       --ref GRCh38.fa \
#       --outdir results/variants \
#       --intervals intervals.list \
#       [--batch-size 50] [--threads 4]
###############################################################################

# Defaults
MODE=""
BAM=""
REF=""
SAMPLE=""
OUTDIR=""
INTERVALS=""
GVCF_DIR=""
BATCH_SIZE=50
THREADS=4
EXTRA_HC=""
EXTRA_GT=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --mode <single|joint> [options]

Modes:
  single    Run GATK HaplotypeCaller on a single sample → gVCF
  joint     Run GenomicsDBImport + GenotypeGVCFs → cohort VCF

Required (both modes):
  --ref <path>           Reference FASTA (with .fai and .dict)
  --outdir <path>        Output directory

Required (single mode):
  --bam <path>           Input BAM file (coordinate-sorted, indexed)
  --sample <name>        Sample name

Required (joint mode):
  --gvcf-dir <path>      Directory containing per-sample .g.vcf.gz files
  --intervals <path>     Genomic intervals file (BED or .list)

Optional:
  --intervals <path>     Restrict calling to these regions (BED for WES targets)
  --batch-size <int>     GenomicsDBImport batch size (default: $BATCH_SIZE)
  --threads <int>        Number of threads (default: $THREADS)
  --extra-hc <args>      Extra HaplotypeCaller arguments (quoted)
  --extra-gt <args>      Extra GenotypeGVCFs arguments (quoted)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)        shift; MODE="$1" ;;
        --bam)         shift; BAM="$1" ;;
        --ref)         shift; REF="$1" ;;
        --sample)      shift; SAMPLE="$1" ;;
        --outdir)      shift; OUTDIR="$1" ;;
        --intervals)   shift; INTERVALS="$1" ;;
        --gvcf-dir)    shift; GVCF_DIR="$1" ;;
        --batch-size)  shift; BATCH_SIZE="$1" ;;
        --threads)     shift; THREADS="$1" ;;
        --extra-hc)    shift; EXTRA_HC="$1" ;;
        --extra-gt)    shift; EXTRA_GT="$1" ;;
        -h|--help)     usage ;;
        *)             echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Common validation
[[ -n "$MODE" ]]   || { echo "ERROR: --mode is required (single or joint)"; usage; }
[[ -n "$REF" ]]    || { echo "ERROR: --ref is required"; usage; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; }
[[ -f "$REF" ]]    || { echo "ERROR: Reference not found: $REF"; exit 1; }

command -v gatk >/dev/null || { echo "ERROR: gatk not found in PATH"; exit 1; }

# Reference checks
[[ -f "${REF}.fai" ]]          || { echo "ERROR: Missing .fai index for reference"; exit 1; }
DICT="${REF%.fa}.dict"
[[ -f "$DICT" ]] || DICT="${REF%.fasta}.dict"
[[ -f "$DICT" ]]               || { echo "ERROR: Missing .dict for reference"; exit 1; }


###############################################################################
# MODE: SINGLE — HaplotypeCaller per sample
###############################################################################
if [[ "$MODE" == "single" ]]; then

    [[ -n "$BAM" ]]    || { echo "ERROR: --bam is required for single mode"; usage; }
    [[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required for single mode"; usage; }
    [[ -f "$BAM" ]]    || { echo "ERROR: BAM not found: $BAM"; exit 1; }

    # Check BAM index
    if [[ ! -f "${BAM}.bai" ]] && [[ ! -f "${BAM%.*}.bai" ]]; then
        echo "ERROR: BAM index not found. Run: samtools index $BAM"
        exit 1
    fi

    GVCF_DIR_OUT="${OUTDIR}/gvcf"
    LOG_DIR="${OUTDIR}/logs"
    mkdir -p "$GVCF_DIR_OUT" "$LOG_DIR"

    GVCF_OUT="${GVCF_DIR_OUT}/${SAMPLE}.g.vcf.gz"
    LOG="${LOG_DIR}/${SAMPLE}.haplotypecaller.log"
    CHECK="${GVCF_DIR_OUT}/${SAMPLE}.hc.done"

    # Build intervals flag
    INTERVAL_FLAG=""
    if [[ -n "$INTERVALS" ]]; then
        [[ -f "$INTERVALS" ]] || { echo "ERROR: Intervals file not found: $INTERVALS"; exit 1; }
        INTERVAL_FLAG="-L $INTERVALS"
    fi

    echo "========================================="
    echo "GATK HaplotypeCaller"
    echo "========================================="
    echo "  Sample:     $SAMPLE"
    echo "  BAM:        $(basename $BAM)"
    echo "  Reference:  $(basename $REF)"
    echo "  Intervals:  ${INTERVALS:-whole genome}"
    echo "  Threads:    $THREADS"
    echo "  Output:     $GVCF_OUT"
    echo "========================================="

    if [[ -f "$CHECK" ]]; then
        echo "HaplotypeCaller already completed — skipping"
    else
        echo "Running HaplotypeCaller..."

        gatk HaplotypeCaller \
            -R "$REF" \
            -I "$BAM" \
            -O "$GVCF_OUT" \
            -ERC GVCF \
            $INTERVAL_FLAG \
            --native-pair-hmm-threads "$THREADS" \
            $EXTRA_HC \
            2> "$LOG"

        # Verify output
        [[ -f "$GVCF_OUT" ]]     || { echo "ERROR: gVCF not created"; exit 1; }
        [[ -f "${GVCF_OUT}.tbi" ]] || { echo "ERROR: gVCF index not created"; exit 1; }

        touch "$CHECK"
        echo "HaplotypeCaller complete: $GVCF_OUT"
    fi

    echo ""
    echo "========================================="
    echo "SAMPLE $SAMPLE — gVCF COMPLETE"
    echo "========================================="
    echo "  gVCF:  $GVCF_OUT"
    echo "  Index: ${GVCF_OUT}.tbi"
    echo "  Log:   $LOG"
    echo "========================================="


###############################################################################
# MODE: JOINT — GenomicsDBImport + GenotypeGVCFs
###############################################################################
elif [[ "$MODE" == "joint" ]]; then

    [[ -n "$GVCF_DIR" ]]  || { echo "ERROR: --gvcf-dir is required for joint mode"; usage; }
    [[ -d "$GVCF_DIR" ]]  || { echo "ERROR: gVCF directory not found: $GVCF_DIR"; exit 1; }
    [[ -n "$INTERVALS" ]] || { echo "ERROR: --intervals is required for joint mode"; usage; }
    [[ -f "$INTERVALS" ]] || { echo "ERROR: Intervals file not found: $INTERVALS"; exit 1; }

    GENOMICSDB="${OUTDIR}/genomicsdb"
    JOINT_DIR="${OUTDIR}/joint"
    LOG_DIR="${OUTDIR}/logs"
    mkdir -p "$JOINT_DIR" "$LOG_DIR"

    COHORT_VCF="${JOINT_DIR}/cohort.vcf.gz"
    LOG_IMPORT="${LOG_DIR}/genomicsdb_import.log"
    LOG_GENOTYPE="${LOG_DIR}/genotype_gvcfs.log"
    CHECK_IMPORT="${OUTDIR}/.genomicsdb_import.done"
    CHECK_GENOTYPE="${OUTDIR}/.genotype_gvcfs.done"

    # Collect all gVCFs
    GVCF_FILES=( "${GVCF_DIR}"/*.g.vcf.gz )
    if [[ ${#GVCF_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No .g.vcf.gz files found in $GVCF_DIR"
        exit 1
    fi

    # Build -V arguments
    GVCF_ARGS=""
    for gvcf in "${GVCF_FILES[@]}"; do
        GVCF_ARGS+=" -V $gvcf"
    done

    echo "========================================="
    echo "GATK Joint Genotyping"
    echo "========================================="
    echo "  gVCFs:      ${#GVCF_FILES[@]} samples"
    echo "  Reference:  $(basename $REF)"
    echo "  Intervals:  $INTERVALS"
    echo "  Batch size: $BATCH_SIZE"
    echo "  Output:     $COHORT_VCF"
    echo "========================================="

    # --- GenomicsDBImport ---
    if [[ -f "$CHECK_IMPORT" ]]; then
        echo "GenomicsDBImport already completed — skipping"
    else
        echo ""
        echo "Step 1/2: GenomicsDBImport..."
        echo "  Importing ${#GVCF_FILES[@]} gVCFs into GenomicsDB workspace"

        # Remove existing workspace if present (GenomicsDBImport requires clean dir)
        [[ -d "$GENOMICSDB" ]] && rm -rf "$GENOMICSDB"

        gatk GenomicsDBImport \
            $GVCF_ARGS \
            --genomicsdb-workspace-path "$GENOMICSDB" \
            -L "$INTERVALS" \
            --batch-size "$BATCH_SIZE" \
            --reader-threads "$THREADS" \
            2> "$LOG_IMPORT"

        [[ -d "$GENOMICSDB" ]] || { echo "ERROR: GenomicsDB workspace not created"; exit 1; }

        touch "$CHECK_IMPORT"
        echo "  GenomicsDBImport complete: $GENOMICSDB"
    fi

    # --- GenotypeGVCFs ---
    if [[ -f "$CHECK_GENOTYPE" ]]; then
        echo "GenotypeGVCFs already completed — skipping"
    else
        echo ""
        echo "Step 2/2: GenotypeGVCFs..."

        gatk GenotypeGVCFs \
            -R "$REF" \
            -V "gendb://$GENOMICSDB" \
            -O "$COHORT_VCF" \
            $EXTRA_GT \
            2> "$LOG_GENOTYPE"

        [[ -f "$COHORT_VCF" ]] || { echo "ERROR: Cohort VCF not created"; exit 1; }

        # Index the output
        if [[ ! -f "${COHORT_VCF}.tbi" ]]; then
            echo "  Indexing cohort VCF..."
            gatk IndexFeatureFile -I "$COHORT_VCF"
        fi

        touch "$CHECK_GENOTYPE"
        echo "  GenotypeGVCFs complete: $COHORT_VCF"
    fi

    echo ""
    echo "========================================="
    echo "JOINT GENOTYPING COMPLETE"
    echo "========================================="
    echo "  Samples:    ${#GVCF_FILES[@]}"
    echo "  Cohort VCF: $COHORT_VCF"
    echo "  Index:      ${COHORT_VCF}.tbi"
    echo "  GenomicsDB: $GENOMICSDB"
    echo ""
    echo "  Logs:"
    echo "    Import:   $LOG_IMPORT"
    echo "    Genotype: $LOG_GENOTYPE"
    echo "========================================="

else
    echo "ERROR: Unknown mode '$MODE'. Use 'single' or 'joint'."
    usage
fi
