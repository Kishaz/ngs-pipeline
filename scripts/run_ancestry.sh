#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# Ancestry Estimation via Singularity Container
# Stages BAM files, handles sample name sanitization, runs batch ancestry
#
# Usage:
#   bash run_ancestry.sh \
#       --bam-dir results/alignment/bam \
#       --outdir results/ancestry \
#       --container /path/to/ancestry-pipeline.sif \
#       --seq-type exome \
#       [--threads 4] [--subpops]
#
# Or single-sample mode:
#   bash run_ancestry.sh \
#       --bam results/alignment/bam/Sample.sorted.markdup.bam \
#       --sample SampleName \
#       --outdir results/ancestry \
#       --container /path/to/ancestry-pipeline.sif \
#       --seq-type wgs
###############################################################################

# Defaults
THREADS=4
SEQ_TYPE="exome"
SUBPOPS=""
BAM=""
BAM_DIR=""
SAMPLE=""
OUTDIR=""
CONTAINER=""
RUNTIME="singularity"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required:
  --outdir <path>                Output directory for ancestry results
  --container <path>             Path to ancestry-pipeline Singularity (.sif) image

Input (one required):
  --bam <path>                   Single BAM file (requires --sample)
  --bam-dir <path>               Directory of BAM files (batch mode)

Optional:
  --sample <name>                Sample name (single-sample mode only)
  --seq-type <wgs|exome|rna>     Sequencing type (default: $SEQ_TYPE)
  --threads <int>                Number of threads (default: $THREADS)
  --subpops                      Run sub-population analysis (K=23)
  --runtime <singularity|docker> Container runtime (default: $RUNTIME)

Notes:
  - BAMs must be aligned to GRCh38 and indexed (.bai)
  - Sample names with underscores are automatically converted to hyphens
    (ADMIXTURE requirement)
  - In batch mode, all BAMs in the directory are processed together
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)        shift; BAM="$1" ;;
        --bam-dir)    shift; BAM_DIR="$1" ;;
        --sample)     shift; SAMPLE="$1" ;;
        --outdir)     shift; OUTDIR="$1" ;;
        --container)  shift; CONTAINER="$1" ;;
        --seq-type)   shift; SEQ_TYPE="$1" ;;
        --threads)    shift; THREADS="$1" ;;
        --subpops)    SUBPOPS="-subpops" ;;
        --runtime)    shift; RUNTIME="$1" ;;
        -h|--help)    usage ;;
        *)            echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Validation
[[ -n "$OUTDIR" ]]    || { echo "ERROR: --outdir is required"; usage; }
[[ -n "$CONTAINER" ]] || { echo "ERROR: --container is required"; usage; }
[[ -f "$CONTAINER" ]] || { echo "ERROR: Container not found: $CONTAINER"; exit 1; }

if [[ -z "$BAM" ]] && [[ -z "$BAM_DIR" ]]; then
    echo "ERROR: --bam or --bam-dir is required"
    usage
fi

if [[ -n "$BAM" ]] && [[ -z "$SAMPLE" ]]; then
    echo "ERROR: --sample is required for single-BAM mode"
    usage
fi

mkdir -p "$OUTDIR"
OUTDIR_ABS=$(cd "$OUTDIR" && pwd -P)

###############################################################################
# STAGING — Symlink BAMs with sanitized names
###############################################################################
STAGING_DIR="${OUTDIR_ABS}/.staging"
mkdir -p "$STAGING_DIR"

stage_bam() {
    local bam_path="$1"
    local sample_name="$2"

    # Sanitize: replace underscores with hyphens (ADMIXTURE requirement)
    local safe_name="${sample_name//_/-}"

    local bam_abs=$(cd "$(dirname "$bam_path")" && pwd -P)/$(basename "$bam_path")
    local bai_abs=""

    # Find index
    if [[ -f "${bam_abs}.bai" ]]; then
        bai_abs="${bam_abs}.bai"
    elif [[ -f "${bam_abs%.*}.bai" ]]; then
        bai_abs="${bam_abs%.*}.bai"
    else
        echo "WARNING: No .bai index for $(basename "$bam_path") — GATK will likely fail"
    fi

    # Symlink with safe name
    local dst_bam="${STAGING_DIR}/${safe_name}.bam"
    local dst_bai="${STAGING_DIR}/${safe_name}.bam.bai"

    [[ -L "$dst_bam" ]] && rm -f "$dst_bam"
    [[ -L "$dst_bai" ]] && rm -f "$dst_bai"

    ln -sf "$bam_abs" "$dst_bam"
    [[ -n "$bai_abs" ]] && ln -sf "$bai_abs" "$dst_bai"

    if [[ "$sample_name" != "$safe_name" ]]; then
        echo "  Staged: ${sample_name} -> ${safe_name} (underscores removed)"
    else
        echo "  Staged: ${sample_name}"
    fi
}

echo "========================================="
echo "Ancestry Estimation Pipeline"
echo "========================================="
echo "  Container: $(basename $CONTAINER)"
echo "  Seq Type:  $SEQ_TYPE"
echo "  Threads:   $THREADS"
echo "  Subpops:   $([ -n "$SUBPOPS" ] && echo 'yes' || echo 'no')"
echo "  Output:    $OUTDIR_ABS"
echo "-----------------------------------------"
echo "Staging BAM files..."

NUM_STAGED=0

if [[ -n "$BAM_DIR" ]]; then
    # Batch mode — stage all BAMs in directory
    BAM_DIR_ABS=$(cd "$BAM_DIR" && pwd -P)

    # If the bam-dir IS the staging directory (pre-staged by Snakemake),
    # skip re-staging to avoid creating duplicate symlinks.
    if [[ "$BAM_DIR_ABS" == "$STAGING_DIR" ]]; then
        echo "  Using pre-staged directory (skipping re-staging)"
        for bamfile in "${BAM_DIR_ABS}"/*.bam; do
            [[ -f "$bamfile" ]] || continue
            echo "  Staged: $(basename "$bamfile" .bam)"
            NUM_STAGED=$((NUM_STAGED + 1))
        done
    else
        for bamfile in "${BAM_DIR_ABS}"/*.bam; do
            [[ -f "$bamfile" ]] || continue
            bname=$(basename "$bamfile" .bam)
            # Strip common suffixes for clean sample names
            bname=$(echo "$bname" | sed 's/\.sorted\.markdup\.recal$//' \
                                  | sed 's/\.sorted\.markdup$//' \
                                  | sed 's/\.sorted$//' \
                                  | sed 's/\.markdup$//')
            stage_bam "$bamfile" "$bname"
            NUM_STAGED=$((NUM_STAGED + 1))
        done
    fi

    INPUT_BIND="${STAGING_DIR}:/input_data:ro"
    PIPELINE_ARGS="-bam_dir /input_data"
else
    # Single-sample mode
    [[ -f "$BAM" ]] || { echo "ERROR: BAM file not found: $BAM"; exit 1; }
    stage_bam "$BAM" "$SAMPLE"
    NUM_STAGED=1

    SAFE_SAMPLE="${SAMPLE//_/-}"
    INPUT_BIND="${STAGING_DIR}:/input_data:ro"
    PIPELINE_ARGS="-bam /input_data/${SAFE_SAMPLE}.bam -sample ${SAFE_SAMPLE}"
fi

echo ""
echo "Staged $NUM_STAGED sample(s)"
echo "-----------------------------------------"

PIPELINE_ARGS="${PIPELINE_ARGS} -outdir /output -seq_type ${SEQ_TYPE} -threads ${THREADS} ${SUBPOPS}"

###############################################################################
# RUN CONTAINER
###############################################################################
echo "Running ancestry pipeline..."
echo "  Runtime: $RUNTIME"
echo ""

if [[ "$RUNTIME" == "singularity" ]]; then
    singularity run \
        --bind "${INPUT_BIND}" \
        --bind "${OUTDIR_ABS}:/output" \
        "$CONTAINER" \
        $PIPELINE_ARGS

elif [[ "$RUNTIME" == "docker" ]]; then
    docker run --rm \
        -v "${INPUT_BIND}" \
        -v "${OUTDIR_ABS}:/output" \
        "$CONTAINER" \
        $PIPELINE_ARGS

else
    echo "ERROR: Unknown runtime '$RUNTIME'. Use 'singularity' or 'docker'."
    exit 1
fi

###############################################################################
# CLEANUP staging
###############################################################################
echo ""
echo "Cleaning up staging directory..."
rm -rf "$STAGING_DIR"

###############################################################################
# VERIFY OUTPUTS
###############################################################################
echo ""
echo "========================================="
echo "ANCESTRY ESTIMATION COMPLETE"
echo "========================================="

if [[ -f "${OUTDIR_ABS}/ancestry_summary_superpops.tsv" ]]; then
    N_SAMPLES=$(tail -n +2 "${OUTDIR_ABS}/ancestry_summary_superpops.tsv" | wc -l | tr -d ' ')
    echo "  Samples processed: $N_SAMPLES"
    echo "  Summary:  ${OUTDIR_ABS}/ancestry_summary_superpops.tsv"
else
    echo "  WARNING: ancestry_summary_superpops.tsv not found"
fi

[[ -f "${OUTDIR_ABS}/ancestry_pca.pdf" ]] && \
    echo "  PCA Plot: ${OUTDIR_ABS}/ancestry_pca.pdf"

if [[ -n "$SUBPOPS" ]] && [[ -f "${OUTDIR_ABS}/ancestry_summary_subpops.tsv" ]]; then
    echo "  Subpops:  ${OUTDIR_ABS}/ancestry_summary_subpops.tsv"
fi

# List per-sample ancestry files
echo ""
echo "  Per-sample results:"
for f in "${OUTDIR_ABS}"/*.ancestry_superpops.tsv; do
    [[ -f "$f" ]] || continue
    sample=$(basename "$f" .ancestry_superpops.tsv)
    echo "    $sample"
done

echo ""
echo "========================================="
