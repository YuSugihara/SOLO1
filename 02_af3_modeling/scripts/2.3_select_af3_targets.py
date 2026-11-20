#!/usr/bin/env python3
"""Step 2.3 script: drop the representatives the AlphaFold Server cannot accept.

Each representative is modelled as a hexamer together with 50 copies of a
ligand, so the size of a job grows with the length of the region. The longest
representatives exceed the size limit of the AlphaFold Server and cannot be
submitted. They are identified here by a length threshold; every
representative appears in the report with its length and the resulting
decision.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-fasta", required=True, type=Path,
                   help="Non-redundant modelled regions.")
    p.add_argument("--output-fasta", required=True, type=Path,
                   help="Regions submitted to AlphaFold 3.")
    p.add_argument("--report-tsv", required=True, type=Path,
                   help="Per-representative report of length and decision.")
    p.add_argument("--max-length", type=int, default=700,
                   help="Longest region that can be submitted, in residues "
                        "(default: %(default)s). Any value between the longest "
                        "accepted and the shortest rejected region reproduces "
                        "the selection.")
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

    records = list(read_fasta(args.input_fasta))
    if len({name for name, _ in records}) != len(records):
        sys.exit(f"error: {args.input_fasta} contains duplicate identifiers")

    kept = excluded = 0
    with open(args.output_fasta, "w") as fasta_out, open(args.report_tsv, "w") as report_out:
        report_out.write("modelled_id\tlength\tdecision\n")
        for name, sequence in records:
            if len(sequence) > args.max_length:
                report_out.write(f"{name}\t{len(sequence)}\texceeds_server_limit\n")
                excluded += 1
            else:
                report_out.write(f"{name}\t{len(sequence)}\tsubmitted\n")
                fasta_out.write(f">{name}\n{sequence}\n")
                kept += 1

    print(f"Step 2.3 select_af3_targets completed: "
          f"representatives={len(records)}, submitted={kept}, "
          f"exceeds_server_limit={excluded}, max_length={args.max_length}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
