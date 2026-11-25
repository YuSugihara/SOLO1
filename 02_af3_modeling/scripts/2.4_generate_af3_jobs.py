#!/usr/bin/env python3
"""Step 2.4 script: write AlphaFold Server job files for the selected regions.

Each region is submitted as a homo-oligomer together with a fixed number of
copies of a ligand, in the batch format the AlphaFold Server accepts. Model
seeds are left empty so that the server assigns an independent seed to every
submission; submitting the same files once per replicate therefore yields
independent predictions.
"""

import argparse
import json
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-fasta", required=True, type=Path,
                   help="Regions to submit, one record per job.")
    p.add_argument("--output-dir", required=True, type=Path,
                   help="Directory to write the job files into.")
    p.add_argument("--output-prefix", default="2.4_af3_jobs",
                   help="Basename of the job files (default: %(default)s).")
    p.add_argument("--jobs-per-file", type=int, default=20,
                   help="Jobs per file; the server accepts up to 20 "
                        "(default: %(default)s).")
    p.add_argument("--copies", type=int, default=6,
                   help="Copies of the protein chain per job (default: %(default)s).")
    p.add_argument("--ligand", default="CCD_OLA",
                   help="Ligand added to every job, as a CCD code "
                        "(default: %(default)s). Pass an empty string for none.")
    p.add_argument("--ligand-count", type=int, default=50,
                   help="Copies of the ligand per job (default: %(default)s).")
    p.add_argument("--no-structure-template", action="store_true",
                   help="Turn off the use of structure templates.")
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
    if args.jobs_per_file < 1:
        sys.exit("error: --jobs-per-file must be at least 1")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for name, sequence in read_fasta(args.input_fasta):
        sequences = [{
            "proteinChain": {
                "sequence": sequence,
                "count": args.copies,
                "useStructureTemplate": not args.no_structure_template,
            }
        }]
        if args.ligand:
            sequences.append({"ligand": {"ligand": args.ligand,
                                         "count": args.ligand_count}})
        jobs.append({
            "name": name,
            "modelSeeds": [],
            "sequences": sequences,
            "dialect": "alphafoldserver",
            "version": 1,
        })

    if not jobs:
        sys.exit(f"error: no records in {args.input_fasta}")

    files = 0
    for start in range(0, len(jobs), args.jobs_per_file):
        files += 1
        path = args.output_dir / f"{args.output_prefix}_{files:03d}.json"
        with open(path, "w") as handle:
            json.dump(jobs[start:start + args.jobs_per_file], handle, indent=2)
            handle.write("\n")

    print(f"Step 2.4 generate_af3_jobs completed: jobs={len(jobs)}, files={files}, "
          f"copies={args.copies}, ligand={args.ligand or 'none'}x{args.ligand_count}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
