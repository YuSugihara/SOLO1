#!/usr/bin/env python3
"""Step 7.1 script: pull the SOLO1 clade out of the full-length NLR set.

The members of the clade were read off the plant-wide NLR tree and are listed
by identifier. This step extracts their full-length sequences, and fails if any
listed identifier is absent, so that a silently truncated clade cannot pass
unnoticed.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--id-list", required=True, type=Path,
                   help="Identifiers of the clade members, one per line.")
    p.add_argument("--input-fasta", required=True, type=Path,
                   help="Full-length NLR sequences to draw from.")
    p.add_argument("--output-fasta", required=True, type=Path,
                   help="Sequences of the listed identifiers.")
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
    wanted = [line.strip() for line in args.id_list.read_text().splitlines() if line.strip()]
    if len(set(wanted)) != len(wanted):
        sys.exit(f"error: {args.id_list} lists the same identifier more than once")

    sequences = {name: seq for name, seq in read_fasta(args.input_fasta)
                 if name in set(wanted)}
    missing = [name for name in wanted if name not in sequences]
    if missing:
        sys.exit(f"error: {len(missing)} listed identifiers are absent from "
                 f"{args.input_fasta}, for example {missing[:3]}")

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_fasta, "w") as out:
        for name in wanted:
            out.write(f">{name}\n{sequences[name]}\n")
    print(f"Step 7.1 extract_clade_sequences completed: sequences={len(wanted)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
