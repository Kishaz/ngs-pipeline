#!/usr/bin/env python3
"""
Generate RNA-seq QC reference files from a GENCODE/Ensembl GTF, WITHOUT
UCSC tools (gtfToGenePred/genePredToBed) or any network access.

Outputs (any subset, depending on which --out flags are given):
  --bed12   BED12 gene model      -> RSeQC infer_experiment.py / read_distribution.py
  --refflat refFlat (genePred+)   -> Picard CollectRnaSeqMetrics REF_FLAT
  --rrna    rRNA interval_list     -> Picard CollectRnaSeqMetrics RIBOSOMAL_INTERVALS

The rRNA interval_list needs the reference sequence dictionary (@SQ lines) so
its contigs match the aligned BAMs; pass it with --dict (the genome .dict, or a
.fai — @SQ lines are synthesized from a .fai if a .dict is not given).

Usage:
  python gtf_to_rnaseq_refs.py \
      --gtf annotation/gencode.v48.primary_assembly.annotation.gtf \
      --dict genome/GRCh38.primary_assembly.genome.dict \
      --bed12 annotation/gencode.v48.genes.bed12 \
      --refflat annotation/gencode.v48.refFlat.txt \
      --rrna annotation/gencode.v48.rRNA.interval_list
"""
import argparse
import re
import sys

# rRNA-like biotypes to mark as ribosomal for Picard
RRNA_TYPES = {"rRNA", "Mt_rRNA", "rRNA_pseudogene"}

_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def parse_attrs(field):
    return dict(_ATTR_RE.findall(field))


def load_gtf(path):
    """
    One streaming pass. Returns:
      tx: transcript_id -> dict(chrom, strand, gene_name, exons=[(s,e)],
                                cds_min, cds_max)   (1-based inclusive coords)
      rrna_genes: list of (chrom, start, end, strand, gene_id)
    """
    tx = {}
    rrna_genes = []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, feat, start, end, _, strand = f[0], f[1], f[2], f[3], f[4], f[5], f[6]
            if feat not in ("exon", "CDS", "gene"):
                continue
            start, end = int(start), int(end)
            attrs = parse_attrs(f[8])
            if feat == "gene":
                gtype = attrs.get("gene_type") or attrs.get("gene_biotype", "")
                if gtype in RRNA_TYPES:
                    rrna_genes.append((chrom, start, end, strand,
                                       attrs.get("gene_id", "rRNA")))
                continue
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            rec = tx.get(tid)
            if rec is None:
                rec = tx[tid] = {
                    "chrom": chrom, "strand": strand,
                    "gene_name": attrs.get("gene_name") or attrs.get("gene_id", tid),
                    "exons": [], "cds_min": None, "cds_max": None,
                }
            if feat == "exon":
                rec["exons"].append((start, end))
            else:  # CDS
                rec["cds_min"] = start if rec["cds_min"] is None else min(rec["cds_min"], start)
                rec["cds_max"] = end if rec["cds_max"] is None else max(rec["cds_max"], end)
    return tx, rrna_genes


def write_bed12(tx, out):
    n = 0
    with open(out, "w") as w:
        for tid, r in tx.items():
            exons = sorted(r["exons"])
            if not exons:
                continue
            chrom_start = exons[0][0] - 1          # BED is 0-based
            chrom_end = exons[-1][1]
            if r["cds_min"] is not None:
                thick_start, thick_end = r["cds_min"] - 1, r["cds_max"]
            else:                                    # non-coding -> zero-length thick
                thick_start = thick_end = chrom_start
            sizes = ",".join(str(e - s + 1) for s, e in exons) + ","
            starts = ",".join(str((s - 1) - chrom_start) for s, e in exons) + ","
            w.write("\t".join(map(str, [
                r["chrom"], chrom_start, chrom_end, tid, 0, r["strand"],
                thick_start, thick_end, 0, len(exons), sizes, starts,
            ])) + "\n")
            n += 1
    return n


def write_refflat(tx, out):
    n = 0
    with open(out, "w") as w:
        for tid, r in tx.items():
            exons = sorted(r["exons"])
            if not exons:
                continue
            tx_start = exons[0][0] - 1
            tx_end = exons[-1][1]
            if r["cds_min"] is not None:
                cds_start, cds_end = r["cds_min"] - 1, r["cds_max"]
            else:
                cds_start = cds_end = tx_end       # convention: no CDS
            estarts = ",".join(str(s - 1) for s, e in exons) + ","
            eends = ",".join(str(e) for s, e in exons) + ","
            w.write("\t".join(map(str, [
                r["gene_name"], tid, r["chrom"], r["strand"],
                tx_start, tx_end, cds_start, cds_end, len(exons), estarts, eends,
            ])) + "\n")
            n += 1
    return n


def sq_lines_from_dict(path):
    with open(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.startswith("@SQ")]


def sq_lines_from_fai(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            name, length = ln.split("\t")[:2]
            out.append(f"@SQ\tSN:{name}\tLN:{length}")
    return out


def write_rrna(rrna_genes, sq_lines, out):
    with open(out, "w") as w:
        w.write("@HD\tVN:1.6\n")
        for sq in sq_lines:
            w.write(sq + "\n")
        for chrom, start, end, strand, gid in rrna_genes:
            w.write(f"{chrom}\t{start}\t{end}\t{strand}\t{gid}\n")
    return len(rrna_genes)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--dict", help="genome .dict (preferred) or .fai for rRNA @SQ header")
    ap.add_argument("--bed12")
    ap.add_argument("--refflat")
    ap.add_argument("--rrna")
    args = ap.parse_args()

    if not (args.bed12 or args.refflat or args.rrna):
        ap.error("nothing to do: give at least one of --bed12 / --refflat / --rrna")

    sys.stderr.write(f"[refs] parsing {args.gtf} ...\n")
    tx, rrna_genes = load_gtf(args.gtf)
    sys.stderr.write(f"[refs] {len(tx)} transcripts, {len(rrna_genes)} rRNA genes\n")

    if args.bed12:
        sys.stderr.write(f"[refs] BED12   -> {args.bed12} ({write_bed12(tx, args.bed12)} transcripts)\n")
    if args.refflat:
        sys.stderr.write(f"[refs] refFlat -> {args.refflat} ({write_refflat(tx, args.refflat)} transcripts)\n")
    if args.rrna:
        if not args.dict:
            ap.error("--rrna requires --dict (genome .dict or .fai) for the @SQ header")
        sq = (sq_lines_from_dict(args.dict) if args.dict.endswith(".dict")
              else sq_lines_from_fai(args.dict))
        if not sq:
            sys.stderr.write(f"[refs] WARNING: no @SQ lines found in {args.dict}\n")
        sys.stderr.write(f"[refs] rRNA    -> {args.rrna} ({write_rrna(rrna_genes, sq, args.rrna)} intervals)\n")

    sys.stderr.write("[refs] done.\n")


if __name__ == "__main__":
    main()
