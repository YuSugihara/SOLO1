#!/usr/bin/env python3
"""Step 2.1 script: cut each NLR-ID down to the region modelled by AlphaFold 3.

For every NLR-ID, the modelled region runs from the N terminus to the end of
the NB-ARC domain, as annotated by NLRtracker. Sequences containing an
ambiguous residue (X) within that region are dropped, because the AlphaFold
Server rejects them.

Identifiers are rewritten as `<accession>_6x<protein>_1-<end>`, with dots
replaced by underscores, where `6x` records that the region is submitted as a
hexamer. Every input protein appears in the report, whether it was kept or
dropped.
"""

import argparse
import sys
from pathlib import Path

NBARC_ROW = ("REGION", "NBARC")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-fasta", required=True, type=Path,
                   help="Full-length NLR-ID protein sequences from the discovery set.")
    p.add_argument("--input-domains", required=True, type=Path,
                   help="NLRtracker domain table of the same proteins.")
    p.add_argument("--output-fasta", required=True, type=Path,
                   help="Modelled regions, one record per retained protein.")
    p.add_argument("--report-tsv", required=True, type=Path,
                   help="Per-protein report of the modelled region and the decision.")
    p.add_argument("--ambiguous-residue", default="X",
                   help="Residue character that disqualifies a sequence "
                        "(default: %(default)s).")
    p.add_argument("--copies", type=int, default=6,
                   help="Number of copies the region is modelled as; recorded in "
                        "the identifier (default: %(default)s).")
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


def read_nbarc_ends(path):
    """Return {`<accession>|<protein>`: end} for the NB-ARC region of each protein."""
    ends = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for required in ("assembly_accession", "seqname", "type", "description", "end"):
            if required not in header:
                fail(f"{path} has no '{required}' column")
        i = {name: header.index(name) for name in
             ("assembly_accession", "seqname", "type", "description", "end")}
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(i.values()):
                fail(f"{path} line {line_number} has too few columns")
            if (fields[i["type"]], fields[i["description"]]) != NBARC_ROW:
                continue
            key = f'{fields[i["assembly_accession"]]}|{fields[i["seqname"]]}'
            try:
                end = int(fields[i["end"]])
            except ValueError:
                fail(f"{path} line {line_number} has a non-integer NB-ARC end")
            if key in ends:
                fail(f"{key} has more than one NB-ARC region; the modelled region "
                     f"would be ambiguous")
            ends[key] = end
    return ends


def main():
    args = parse_args()
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.report_tsv.parent.mkdir(parents=True, exist_ok=True)

    ends = read_nbarc_ends(args.input_domains)
    kept = dropped = 0

    with open(args.output_fasta, "w") as fasta_out, open(args.report_tsv, "w") as report_out:
        report_out.write("assembly_accession\tprotein_id\tmodelled_id\t"
                         "full_length\tnbarc_end\tdecision\treason\n")
        for name, sequence in read_fasta(args.input_fasta):
            accession, _, protein = name.partition("|")
            if name not in ends:
                fail(f"{name}: no NB-ARC region in {args.input_domains}")
            end = ends[name]
            if end > len(sequence):
                fail(f"{name}: NB-ARC ends at {end} but the sequence is "
                     f"{len(sequence)} residues long")
            region = sequence[:end]
            modelled_id = (f'{accession.replace(".", "_")}_{args.copies}x'
                           f'{protein.replace(".", "_")}_1-{end}')
            if args.ambiguous_residue in region:
                decision, reason = "dropped", f"ambiguous_residue_{args.ambiguous_residue}"
                dropped += 1
            else:
                decision, reason = "kept", "modelled"
                kept += 1
                fasta_out.write(f">{modelled_id}\n{region}\n")
            report_out.write(f"{accession}\t{protein}\t{modelled_id}\t{len(sequence)}\t"
                             f"{end}\t{decision}\t{reason}\n")

    print(f"Step 2.1 extract_modelled_region completed: "
          f"kept={kept}, dropped={dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
