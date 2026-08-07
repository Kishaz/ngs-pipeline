#!/usr/bin/env python3
"""
Auto-generate samples.tsv from a directory of FASTQ or BAM files.

Scans input directory for .fastq.gz and .bam files, detects PE/SE from
standard Illumina naming patterns, auto-detects BAM processing stage,
and writes a pipeline-compatible samples.tsv.

Usage (standalone):
    python generate_samples.py \
        --input-dir /path/to/data \
        --seq-type rnaseq \
        --output config/samples.tsv

Usage (from Snakefile — called automatically when input_dir is set in config).
"""

import argparse
import glob
import os
import re
import sys

# ── FASTQ R1 patterns (order = priority) ─────────────────────────────────────
# Each tuple: (R1 regex, R1→R2 substitution pairs)
R1_PATTERNS = [
    (r"_R1_001\.fastq\.gz$", ("_R1_001.", "_R2_001.")),
    (r"_R1\.fastq\.gz$",     ("_R1.",     "_R2.")),
    (r"_1\.fastq\.gz$",      ("_1.",      "_2.")),
    (r"\.R1\.fastq\.gz$",    (".R1.",     ".R2.")),
]

# ── BAM stage suffixes (longest match first) ──────────────────────────────────
BAM_STAGES = [
    (".sorted.markdup.recal.bam", "recal"),
    (".sorted.markdup.bam",       "markdup"),
    (".sorted.bam",               "sorted"),
    (".raw.bam",                  "raw"),
]


def discover_fastqs(input_dir):
    """Discover FASTQ samples, detect PE/SE from R1/R2 patterns."""
    files = sorted(glob.glob(os.path.join(input_dir, "*.fastq.gz")))
    if not files:
        return []

    samples = []
    seen_r2 = set()

    for fpath in files:
        fname = os.path.basename(fpath)
        abs_path = os.path.abspath(fpath)

        if abs_path in seen_r2:
            continue

        for pattern, (r1_str, r2_str) in R1_PATTERNS:
            if re.search(pattern, fname):
                sample_name = re.sub(pattern, "", fname)
                r2_path = abs_path.replace(r1_str, r2_str)
                if os.path.isfile(r2_path):
                    samples.append({
                        "sample": sample_name,
                        "R1": abs_path,
                        "R2": r2_path,
                        "bam": "",
                        "bam_stage": "",
                    })
                    seen_r2.add(r2_path)
                else:
                    samples.append({
                        "sample": sample_name,
                        "R1": abs_path,
                        "R2": "",
                        "bam": "",
                        "bam_stage": "",
                    })
                break

    return samples


def discover_bams(input_dir):
    """Discover BAM samples, auto-detect processing stage from filename."""
    files = sorted(glob.glob(os.path.join(input_dir, "*.bam")))
    if not files:
        return []

    samples = []
    for fpath in files:
        fname = os.path.basename(fpath)
        abs_path = os.path.abspath(fpath)

        stage = "raw"
        sample_name = fname.replace(".bam", "")

        for suffix, stage_name in BAM_STAGES:
            if fname.endswith(suffix):
                sample_name = fname[: -len(suffix)]
                stage = stage_name
                break

        # Check for index
        has_index = os.path.isfile(abs_path + ".bai") or os.path.isfile(
            abs_path.replace(".bam", ".bai")
        )
        if not has_index:
            print(
                f"  WARNING: No .bai index for {fname} — pipeline will create one",
                file=sys.stderr,
            )

        samples.append({
            "sample": sample_name,
            "R1": "",
            "R2": "",
            "bam": abs_path,
            "bam_stage": stage,
        })

    return samples


def validate_samples(samples):
    """Check for duplicates and readability."""
    seen = {}
    errors = []

    for s in samples:
        name = s["sample"]
        source = s.get("bam") or s.get("R1", "")
        if name in seen:
            errors.append(
                f"  Duplicate sample name '{name}':\n"
                f"    1) {seen[name]}\n"
                f"    2) {source}"
            )
        seen[name] = source

        # Check files are readable
        for key in ("R1", "R2", "bam"):
            path = s.get(key, "")
            if path and not os.access(path, os.R_OK):
                errors.append(f"  File not readable: {path}")

    if errors:
        print("ERROR: Validation failed:", file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)


def write_tsv(samples, seq_type, output_path):
    """Write samples.tsv."""
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with open(output_path, "w") as f:
        f.write("sample\tR1\tR2\tseq_type\tbam\tbam_stage\n")
        for s in samples:
            f.write(
                f"{s['sample']}\t"
                f"{s['R1']}\t"
                f"{s['R2']}\t"
                f"{seq_type}\t"
                f"{s['bam']}\t"
                f"{s['bam_stage']}\n"
            )


def print_summary(samples):
    """Print discovery summary to stderr."""
    fastq_samples = [s for s in samples if not s["bam"]]
    bam_samples = [s for s in samples if s["bam"]]
    pe = sum(1 for s in fastq_samples if s["R2"])
    se = len(fastq_samples) - pe

    print("─" * 50, file=sys.stderr)
    print("  Sample Discovery Summary", file=sys.stderr)
    print("─" * 50, file=sys.stderr)
    print(f"  Total samples:   {len(samples)}", file=sys.stderr)

    if fastq_samples:
        print(f"  FASTQ samples:   {len(fastq_samples)}", file=sys.stderr)
        print(f"    Paired-end:    {pe}", file=sys.stderr)
        print(f"    Single-end:    {se}", file=sys.stderr)

    if bam_samples:
        print(f"  BAM samples:     {len(bam_samples)}", file=sys.stderr)
        stages = {}
        for s in bam_samples:
            stage = s["bam_stage"]
            stages[stage] = stages.get(stage, 0) + 1
        for stage, count in sorted(stages.items()):
            print(f"    {stage:15s} {count}", file=sys.stderr)

    print("─" * 50, file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Auto-generate samples.tsv from FASTQ/BAM directory"
    )
    parser.add_argument(
        "--input-dir", required=True, help="Directory containing FASTQ or BAM files"
    )
    parser.add_argument(
        "--seq-type",
        required=True,
        choices=["rnaseq", "wes", "wgs"],
        help="Sequencing type for all samples",
    )
    parser.add_argument(
        "--output",
        default="config/samples.tsv",
        help="Output TSV path (default: config/samples.tsv)",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    # Discover samples
    fastq_samples = discover_fastqs(args.input_dir)
    bam_samples = discover_bams(args.input_dir)
    all_samples = fastq_samples + bam_samples

    if not all_samples:
        print(
            f"ERROR: No .fastq.gz or .bam files found in {args.input_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Validate
    validate_samples(all_samples)

    # Summary
    print_summary(all_samples)

    # Write
    write_tsv(all_samples, args.seq_type, args.output)
    print(f"  Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
