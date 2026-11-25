#!/usr/bin/env python3
"""Step 1.2 script: remove candidates whose extra domain is an annotation artefact.

Two classes of candidate carry an extra domain that does not represent an
integrated domain, and are removed:

  transposon   the protein carries a domain annotated as TRANSPOSON and is
               therefore transposon-derived rather than an immune receptor;
  lrr_internal every extra domain of the protein falls inside the LRR region,
               that is, between the start of the first and the end of the last
               annotated LRR. Such a domain is a repeat miscalled as a separate
               domain, not a domain integrated at a receptor terminus.

A protein is retained when it carries at least one extra domain outside the LRR
region and no TRANSPOSON domain. The report lists every candidate once, with
the decision and the evidence behind it.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

LRR_REGION = ("REGION", "LRR")
DOMAIN = "DOMAIN"
TRANSPOSON = "TRANSPOSON"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Directory holding the Step 1.1 output files.")
    p.add_argument("--output-dir", required=True, type=Path,
                   help="Directory to write the filtered files into.")
    p.add_argument("--report-tsv", required=True, type=Path,
                   help="Per-protein report of the retain or remove decision.")
    p.add_argument("--input-prefix", default="1.1_nlr_id_candidates",
                   help="Basename shared by the Step 1.1 files (default: %(default)s).")
    p.add_argument("--output-prefix", default="1.2_nlr_id_filtered",
                   help="Basename shared by the output files (default: %(default)s).")
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


def read_gff_rows(path):
    """Return {`<accession>|<protein>`: [row, ...]} from the Step 1.1 GFF."""
    rows = defaultdict(list)
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            seqid = line.split("\t", 1)[0]
            if "|" not in seqid:
                fail(f"{path} line {line_number} has no <accession>|<protein> seqid")
            rows[seqid].append(line)
    return rows


def read_domain_rows(path):
    """Return (header, {`<accession>|<protein>`: [row, ...]}) from the Step 1.1 table.

    The first column is the assembly accession and the second is the protein
    identifier; rows are returned without the accession column.
    """
    rows = defaultdict(list)
    with open(path) as handle:
        header = handle.readline().rstrip("\n")
        if not header.startswith("assembly_accession\t"):
            fail(f"{path} does not start with an assembly_accession column")
        for line_number, line in enumerate(handle, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t", 2)
            if len(fields) < 3:
                fail(f"{path} line {line_number} has too few columns")
            rows[f"{fields[0]}|{fields[1]}"].append(line.split("\t", 1)[1])
    return header, rows


def domain_features(rows, identifier):
    """Return (lrr_span, extra_domains) parsed from the NLRtracker domain rows."""
    lrr_starts, lrr_ends, extra = [], [], []
    for row in rows:
        fields = row.split("\t")
        if len(fields) < 5:
            fail(f"{identifier}: domain row has fewer than 5 columns: {row!r}")
        kind, name = fields[1], fields[2]
        try:
            start, end = int(fields[3]), int(fields[4])
        except ValueError:
            fail(f"{identifier}: domain row has non-integer coordinates: {row!r}")
        if (kind, name) == LRR_REGION:
            lrr_starts.append(start)
            lrr_ends.append(end)
        elif kind == DOMAIN:
            extra.append((name, start, end))
    if not lrr_starts:
        fail(f"{identifier}: no LRR REGION rows; cannot decide whether its extra "
             f"domains lie inside the LRR")
    if not extra:
        fail(f"{identifier}: no DOMAIN rows, yet it was selected as carrying an "
             f"extra domain")
    return (min(lrr_starts), max(lrr_ends)), sorted(extra, key=lambda d: (d[1], d[2], d[0]))


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_tsv.parent.mkdir(parents=True, exist_ok=True)

    src, dst = args.input_dir, args.output_dir
    ip, op = args.input_prefix, args.output_prefix

    proteins = list(read_fasta(src / f"{ip}.fasta"))
    if len({name for name, _ in proteins}) != len(proteins):
        fail(f"{src / (ip + '.fasta')} contains duplicate identifiers")

    nbarc = defaultdict(list)
    for name, seq in read_fasta(src / f"{ip}_nbarc.fasta"):
        accession, _, record = name.partition("|")
        nbarc[f"{accession}|{record.split('_(')[0]}"].append((name, seq))

    gff_rows = read_gff_rows(src / f"{ip}.gff")
    domains_header, domain_rows = read_domain_rows(src / f"{ip}_domains.tsv")

    counts = defaultdict(int)
    with open(dst / f"{op}.fasta", "w") as protein_out, \
         open(dst / f"{op}_nbarc.fasta", "w") as nbarc_out, \
         open(dst / f"{op}.gff", "w") as gff_out, \
         open(dst / f"{op}_domains.tsv", "w") as domains_out, \
         open(args.report_tsv, "w") as report_out:

        domains_out.write(f"{domains_header}\n")
        report_out.write("assembly_accession\tprotein_id\tdecision\treason\t"
                         "lrr_start\tlrr_end\tn_extra_domains\t"
                         "n_extra_domains_inside_lrr\textra_domains\n")

        for name, sequence in proteins:
            accession, _, protein = name.partition("|")
            rows = domain_rows.get(name)
            if not rows:
                fail(f"{name}: no rows in {ip}_domains.tsv")
            (lrr_start, lrr_end), extra = domain_features(rows, name)

            inside = [d for d in extra if lrr_start <= d[1] and d[2] <= lrr_end]
            has_transposon = any(d[0] == TRANSPOSON for d in extra)
            all_inside = len(inside) == len(extra)

            if has_transposon:
                decision, reason = "removed", "transposon"
            elif all_inside:
                decision, reason = "removed", "lrr_internal"
            else:
                decision, reason = "retained", "integrated_domain"
            counts[reason] += 1

            report_out.write(
                f"{accession}\t{protein}\t{decision}\t{reason}\t{lrr_start}\t{lrr_end}\t"
                f"{len(extra)}\t{len(inside)}\t"
                + ";".join(f"{n}:{s}-{e}" for n, s, e in extra) + "\n")

            if decision != "retained":
                continue

            protein_out.write(f">{name}\n{sequence}\n")
            records = nbarc.get(name)
            if not records:
                fail(f"{name}: no NB-ARC sequence")
            for record_name, seq in sorted(records):
                nbarc_out.write(f">{record_name}\n{seq}\n")
            for row in gff_rows.get(name, []):
                gff_out.write(f"{row}\n")
            for row in rows:
                domains_out.write(f"{accession}\t{row}\n")

    retained = counts["integrated_domain"]
    print(f"Step 1.2 filter_nlr_id_candidates completed: "
          f"candidates={len(proteins)}, removed_transposon={counts['transposon']}, "
          f"removed_lrr_internal={counts['lrr_internal']}, retained={retained}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
