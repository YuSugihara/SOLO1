#!/usr/bin/env python3
"""Step 7.3 script: drop sequences too short to align over their full length.

A receptor whose annotation is truncated contributes mostly gaps to a
full-length alignment and destabilises the tree, so sequences below a minimum
length are set aside. Every sequence is recorded with its length and the
decision.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-fasta", required=True, type=Path)
    p.add_argument("--output-fasta", required=True, type=Path)
    p.add_argument("--report-tsv", required=True, type=Path)
    p.add_argument("--min-length", type=int, default=800,
                   help="Shortest sequence kept, in residues (default: %(default)s).")
    return p.parse_args()


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main():
    args = parse_args()
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.report_tsv.parent.mkdir(parents=True, exist_ok=True)
    kept = dropped = 0
    with open(args.output_fasta, "w") as out, open(args.report_tsv, "w") as report:
        report.write("sequence_id\tlength\tdecision\n")
        for name, sequence in read_fasta(args.input_fasta):
            if len(sequence) >= args.min_length:
                out.write(f">{name}\n{sequence}\n")
                report.write(f"{name}\t{len(sequence)}\tkept\n")
                kept += 1
            else:
                report.write(f"{name}\t{len(sequence)}\ttoo_short\n")
                dropped += 1
    print(f"Step 7.3 filter_by_min_length completed: kept={kept}, "
          f"too_short={dropped}, min_length={args.min_length}", file=sys.stderr)


if __name__ == "__main__":
    main()
