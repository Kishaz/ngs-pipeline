#!/usr/bin/env python3
"""
Parse RSeQC infer_experiment.py output and emit the featureCounts -s value.

featureCounts -s:  0 = unstranded, 1 = stranded (fwd), 2 = reverse-stranded.

RSeQC reports the fraction of reads consistent with each strand convention.
We map:
    fwd = fraction explained by  "1++,1--,2+-,2-+"  (PE)  or  "++,--"  (SE)
    rev = fraction explained by  "1+-,1-+,2++,2--"  (PE)  or  "+-,-+"  (SE)
    p_fwd = fwd / (fwd + rev)
      p_fwd >= threshold      -> 1 (forward-stranded)
      p_fwd <= 1 - threshold  -> 2 (reverse-stranded)
      otherwise               -> 0 (unstranded)

Prints ONLY the integer to stdout (so a rule can `cat` it into featureCounts).
A one-line human-readable rationale goes to stderr.
"""
import argparse
import re
import sys

FWD_RE = re.compile(r'explained by "(?:1\+\+,1--,2\+-,2-\+|\+\+,--)":\s*([0-9.]+)')
REV_RE = re.compile(r'explained by "(?:1\+-,1-\+,2\+\+,2--|\+-,-\+)":\s*([0-9.]+)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", help="RSeQC infer_experiment.py output file")
    ap.add_argument("--threshold", type=float, default=0.80,
                    help="fraction to call stranded (default 0.80)")
    args = ap.parse_args()

    text = open(args.infile).read()
    fwd = FWD_RE.search(text)
    rev = REV_RE.search(text)
    if not fwd or not rev:
        sys.stderr.write("[strand] could not parse infer_experiment output; "
                         "defaulting to unstranded (-s 0)\n")
        print(0)
        return

    fwd, rev = float(fwd.group(1)), float(rev.group(1))
    total = fwd + rev
    p_fwd = fwd / total if total > 0 else 0.5

    if p_fwd >= args.threshold:
        s, label = 1, "forward-stranded"
    elif p_fwd <= (1 - args.threshold):
        s, label = 2, "reverse-stranded"
    else:
        s, label = 0, "unstranded"

    sys.stderr.write(
        f"[strand] fwd={fwd:.3f} rev={rev:.3f} p_fwd={p_fwd:.3f} "
        f"-> {label} (featureCounts -s {s})\n"
    )
    print(s)


if __name__ == "__main__":
    main()
