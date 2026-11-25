#!/usr/bin/env python3
"""Step 6.3 script: merge the NB-ARC sets and keep the well-sized sequences.

The NB-ARC sequences of the selected genomes and of the reference set are
combined into the alignment input. Parentheses are stripped from identifiers,
because several downstream tools treat them as Newick syntax, and sequences
outside the expected length range of an NB-ARC domain are dropped: shorter ones
are partial and longer ones carry flanking sequence, and either would distort
the alignment.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-fasta", required=True, nargs="+", type=Path,
                   help="FASTA files to merge, in the order given.")
    p.add_argument("--output-fasta", required=True, type=Path,
                   help="Merged and length-filtered sequences.")
    p.add_argument("--report-tsv", required=True, type=Path,
                   help="Per-sequence report of length and decision.")
    p.add_argument("--min-length", type=int, default=250,
                   help="Shortest sequence kept (default: %(default)s).")
    p.add_argument("--max-length", type=int, default=400,
                   help="Longest sequence kept (default: %(default)s).")
    return p.parse_args()


def fail(message):
    sys.exit(f"error: {message}")


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
    if args.min_length > args.max_length:
        fail("--min-length is greater than --max-length")
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.report_tsv.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    kept = too_short = too_long = 0
    with open(args.output_fasta, "w") as out, open(args.report_tsv, "w") as report:
        report.write("sequence_id\tsource\tlength\tdecision\n")
        for path in args.input_fasta:
            for name, sequence in read_fasta(path):
                identifier = name.replace("(", "").replace(")", "")
                if identifier in seen:
                    fail(f"{identifier} occurs more than once across the inputs")
                seen.add(identifier)
                length = len(sequence)
                if length < args.min_length:
                    decision = "too_short"
                    too_short += 1
                elif length > args.max_length:
                    decision = "too_long"
                    too_long += 1
                else:
                    decision = "kept"
                    kept += 1
                    out.write(f">{identifier}\n{sequence}\n")
                report.write(f"{identifier}\t{path.name}\t{length}\t{decision}\n")

    print(f"Step 6.3 merge_and_filter_nbarc completed: input={len(seen)}, "
          f"kept={kept}, too_short={too_short}, too_long={too_long}", file=sys.stderr)


if __name__ == "__main__":
    main()
