#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Delete a configurable set of feature columns from the aggregated feature table.

This step used to be done manually. The original pipeline kept 29 features per
variant after extraction, but 6 of them (AF, dp and the four likelihood scores)
were removed by hand before adding the annotation features. This script replaces
that manual step.

Usage:
    python3 delete_features.py --input Total_output_feature_0.1.txt \
                               --output feature.txt \
                               [--drop AF dp mosaic_likelihood het_likelihood refhom_likelihood althom_likelihood]
"""

import argparse
import sys

DEFAULT_DROP = [
    "AF",
    "dp",
    "mosaic_likelihood",
    "het_likelihood",
    "refhom_likelihood",
    "althom_likelihood",
]


def main():
    parser = argparse.ArgumentParser(
        description="Delete feature columns from the aggregated feature table."
    )
    parser.add_argument("--input", required=True, help="aggregated feature file (e.g. Total_output_feature_0.1.txt)")
    parser.add_argument("--output", required=True, help="output feature file (e.g. feature.txt)")
    parser.add_argument(
        "--drop",
        nargs="*",
        default=DEFAULT_DROP,
        help="feature columns to delete (default: %s)" % ", ".join(DEFAULT_DROP),
    )
    args = parser.parse_args()

    drop_set = set(args.drop)
    n_data = 0

    with open(args.input, "r", encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8", newline="\n") as fout:
        header_line = fin.readline()
        if not header_line:
            sys.exit("ERROR: input file %s is empty" % args.input)
        header = header_line.rstrip("\r\n").split("\t")

        missing = sorted(c for c in drop_set if c not in header)
        if missing:
            print("WARNING: requested columns not found in header (skipped): %s" % ", ".join(missing),
                  file=sys.stderr)

        keep_idx = [i for i, col in enumerate(header) if col not in drop_set]
        fout.write("\t".join(header[i] for i in keep_idx) + "\n")

        for line in fin:
            line = line.rstrip("\r\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) != len(header):
                # tolerate ragged rows (pad or trim to header length)
                if len(cols) < len(header):
                    cols = cols + [""] * (len(header) - len(cols))
                else:
                    cols = cols[: len(header)]
            fout.write("\t".join(cols[i] for i in keep_idx) + "\n")
            n_data += 1

    print("Deleted %d column(s) (%s); wrote header + %d data rows to %s" %
          (len(drop_set) - len(missing), ", ".join(sorted(drop_set - set(missing))),
           n_data, args.output))


if __name__ == "__main__":
    main()
