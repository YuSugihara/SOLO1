#!/usr/bin/env python3
"""Step 5.1 script: collect the NB-ARC sequences of the modelled candidates.

The phylogeny is built from the NB-ARC domain alone, which is the part of the
receptor that can be aligned across the whole set. For each candidate submitted
to AlphaFold 3, the NB-ARC sequence annotated by NLRtracker is taken from the
discovery set and relabelled as `<species>|<assembly>|<protein>`, so that tips
read as names rather than as accessions.
"""

import argparse
import re
import sys
from pathlib import Path

TARGET_ID = re.compile(r"^(?P<assembly>GC[AF]_\d+_\d+)_\d+x(?P<protein>[A-Z]{2}_\d+_\d+)_")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target-fasta", required=True, type=Path,
                   help="Regions submitted to AlphaFold 3; defines which "
                        "candidates enter the phylogeny.")
    p.add_argument("--nbarc-fasta", required=True, type=Path,
                   help="NB-ARC sequences of the discovery set.")
    p.add_argument("--assemblies-tsv", required=True, type=Path,
                   help="Assembly accessions with their species names.")
    p.add_argument("--output-fasta", required=True, type=Path,
                   help="Relabelled NB-ARC sequences, ready to align.")
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


def read_species(path):
    species = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for column in ("assembly_accession", "species"):
            if column not in header:
                fail(f"{path} has no '{column}' column")
        i_accession, i_species = header.index("assembly_accession"), header.index("species")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            # Keep the genus and species only, joined by an underscore, so that
            # the label stays readable as a tree tip.
            name = "_".join(fields[i_species].split()[:2])
            species[fields[i_accession]] = name
    return species


def main():
    args = parse_args()
    species = read_species(args.assemblies_tsv)

    nbarc = {}
    for name, sequence in read_fasta(args.nbarc_fasta):
        assembly, _, record = name.partition("|")
        protein = record.split("_(")[0]
        key = (assembly, protein)
        if key in nbarc:
            fail(f"{assembly}|{protein} has more than one NB-ARC sequence; "
                 f"the tip label would be ambiguous")
        nbarc[key] = sequence

    written = 0
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_fasta, "w") as out:
        for name, _ in read_fasta(args.target_fasta):
            match = TARGET_ID.match(name)
            if match is None:
                fail(f"cannot read an assembly and protein accession from {name!r}")
            assembly = match.group("assembly")[::-1].replace("_", ".", 1)[::-1]
            protein = match.group("protein")[::-1].replace("_", ".", 1)[::-1]
            if assembly not in species:
                fail(f"{assembly} is not listed in {args.assemblies_tsv}")
            if (assembly, protein) not in nbarc:
                fail(f"no NB-ARC sequence for {assembly}|{protein}")
            out.write(f">{species[assembly]}|{assembly}|{protein}\n"
                      f"{nbarc[(assembly, protein)]}\n")
            written += 1

    print(f"Step 5.1 prepare_nbarc_alignment_input completed: sequences={written}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
