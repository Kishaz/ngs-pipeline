#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# featureCounts — Gene-level read quantification (RNA-seq)
#
# Modes:
#   --mode single    Count reads for one BAM → per-sample counts file
#   --mode merge     Merge per-sample counts into a gene counts matrix
#
# Usage:
#   # Per-sample counting
#   bash run_featurecounts.sh --mode single \
#       --bam sample.sorted.markdup.bam \
#       --gtf gencode.v48.annotation.gtf \
#       --sample SampleName \
#       --outdir results/counts \
#       [--strandedness 0] [--paired] [--threads 4]
#
#   # Merge into matrix
#   bash run_featurecounts.sh --mode merge \
#       --counts-dir results/counts \
#       --outdir results/counts
###############################################################################

# Defaults
MODE=""
BAM=""
GTF=""
SAMPLE=""
OUTDIR=""
STRANDEDNESS=0
PAIRED=""
THREADS=4
COUNTS_DIR=""
EXTRA="-t exon -g gene_id --primary"

usage() {
    cat <<EOF
Usage: $(basename "$0") --mode <single|merge> [options]

Modes:
  single    Run featureCounts on a single BAM
  merge     Merge per-sample count files into a gene counts matrix

Required (single mode):
  --bam <path>             Input BAM file (coordinate-sorted, indexed)
  --gtf <path>             Gene annotation GTF
  --sample <name>          Sample name
  --outdir <path>          Output directory

Required (merge mode):
  --counts-dir <path>      Directory containing *.featureCounts.txt files
  --outdir <path>          Output directory for merged matrix

Optional:
  --strandedness <0|1|2>   0=unstranded, 1=stranded, 2=reverse (default: $STRANDEDNESS)
  --paired                 Enable paired-end counting (--countReadPairs)
  --threads <int>          Number of threads (default: $THREADS)
  --extra <args>           Extra featureCounts arguments (default: "$EXTRA")
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)           shift; MODE="$1" ;;
        --bam)            shift; BAM="$1" ;;
        --gtf)            shift; GTF="$1" ;;
        --sample)         shift; SAMPLE="$1" ;;
        --outdir)         shift; OUTDIR="$1" ;;
        --strandedness)   shift; STRANDEDNESS="$1" ;;
        --paired)         PAIRED="-p --countReadPairs" ;;
        --threads)        shift; THREADS="$1" ;;
        --counts-dir)     shift; COUNTS_DIR="$1" ;;
        --extra)          shift; EXTRA="$1" ;;
        -h|--help)        usage ;;
        *)                echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Common validation
[[ -n "$MODE" ]]   || { echo "ERROR: --mode is required (single or merge)"; usage; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; }
mkdir -p "$OUTDIR"


###############################################################################
# MODE: SINGLE — featureCounts per sample
###############################################################################
if [[ "$MODE" == "single" ]]; then

    [[ -n "$BAM" ]]    || { echo "ERROR: --bam is required for single mode"; usage; }
    [[ -n "$GTF" ]]    || { echo "ERROR: --gtf is required for single mode"; usage; }
    [[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required for single mode"; usage; }
    [[ -f "$BAM" ]]    || { echo "ERROR: BAM not found: $BAM"; exit 1; }
    [[ -f "$GTF" ]]    || { echo "ERROR: GTF not found: $GTF"; exit 1; }

    command -v featureCounts >/dev/null || { echo "ERROR: featureCounts not found (subread)"; exit 1; }

    COUNTS_OUT="${OUTDIR}/${SAMPLE}.featureCounts.txt"
    SUMMARY="${COUNTS_OUT}.summary"
    LOG="${OUTDIR}/${SAMPLE}.featurecounts.log"

    echo "========================================="
    echo "featureCounts — RNA Quantification"
    echo "========================================="
    echo "  Sample:        $SAMPLE"
    echo "  BAM:           $(basename $BAM)"
    echo "  GTF:           $(basename $GTF)"
    echo "  Strandedness:  $STRANDEDNESS (0=unstranded, 1=stranded, 2=reverse)"
    echo "  Paired-end:    $([ -n "$PAIRED" ] && echo 'yes' || echo 'no')"
    echo "  Threads:       $THREADS"
    echo "  Output:        $COUNTS_OUT"
    echo "========================================="

    featureCounts \
        -a "$GTF" \
        -o "$COUNTS_OUT" \
        $EXTRA \
        -s "$STRANDEDNESS" \
        $PAIRED \
        -T "$THREADS" \
        "$BAM" \
        2> "$LOG"

    [[ -f "$COUNTS_OUT" ]] || { echo "ERROR: Counts file not created"; exit 1; }

    # Quick summary
    ASSIGNED=$(grep "Assigned" "$SUMMARY" | awk '{print $2}')
    UNASSIGNED=$(grep "Unassigned_NoFeatures" "$SUMMARY" | awk '{print $2}')
    echo ""
    echo "  Assigned reads:   ${ASSIGNED:-N/A}"
    echo "  Unassigned (no feature): ${UNASSIGNED:-N/A}"

    echo ""
    echo "========================================="
    echo "SAMPLE $SAMPLE — QUANTIFICATION COMPLETE"
    echo "========================================="
    echo "  Counts:  $COUNTS_OUT"
    echo "  Summary: $SUMMARY"
    echo "  Log:     $LOG"
    echo "========================================="


###############################################################################
# MODE: MERGE — Combine per-sample counts into a matrix
###############################################################################
elif [[ "$MODE" == "merge" ]]; then

    [[ -n "$COUNTS_DIR" ]] || { echo "ERROR: --counts-dir is required for merge mode"; usage; }
    [[ -d "$COUNTS_DIR" ]] || { echo "ERROR: Counts directory not found: $COUNTS_DIR"; exit 1; }

    COUNT_FILES=( "${COUNTS_DIR}"/*.featureCounts.txt )
    if [[ ${#COUNT_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No *.featureCounts.txt files found in $COUNTS_DIR"
        exit 1
    fi

    MATRIX_OUT="${OUTDIR}/gene_counts_matrix.tsv"

    echo "========================================="
    echo "Merge Gene Counts Matrix"
    echo "========================================="
    echo "  Input:   ${#COUNT_FILES[@]} count files"
    echo "  Output:  $MATRIX_OUT"
    echo "========================================="

    python3 - "$COUNTS_DIR" "$MATRIX_OUT" << 'PYEOF'
import sys, os, glob
import csv

counts_dir = sys.argv[1]
matrix_out = sys.argv[2]

count_files = sorted(glob.glob(os.path.join(counts_dir, "*.featureCounts.txt")))
if not count_files:
    print("ERROR: No count files found")
    sys.exit(1)

# Parse each file: skip comment lines (#), header is first non-comment line
# featureCounts columns: Geneid, Chr, Start, End, Strand, Length, <count>
samples = {}
gene_order = []

for cf in count_files:
    sample_name = os.path.basename(cf).replace(".featureCounts.txt", "")
    counts = {}

    with open(cf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[0] == "Geneid":
                continue  # header
            gene_id = parts[0]
            count = int(parts[-1])
            counts[gene_id] = count

            if gene_id not in gene_order:
                gene_order.append(gene_id)

    samples[sample_name] = counts
    print(f"  Parsed: {sample_name} ({len(counts)} genes)")

# Write merged matrix
sample_names = sorted(samples.keys())
with open(matrix_out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Geneid"] + sample_names)
    for gene in gene_order:
        row = [gene] + [samples[s].get(gene, 0) for s in sample_names]
        w.writerow(row)

print(f"\nGene counts matrix: {matrix_out}")
print(f"  Genes:   {len(gene_order)}")
print(f"  Samples: {len(sample_names)}")
PYEOF

    echo ""
    echo "========================================="
    echo "MERGE COMPLETE"
    echo "========================================="
    echo "  Matrix: $MATRIX_OUT"
    echo "========================================="

else
    echo "ERROR: Unknown mode '$MODE'. Use 'single' or 'merge'."
    usage
fi
