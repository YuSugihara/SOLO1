#!/usr/bin/env python3
"""Build iTOL LABELS file from FASTA sequence IDs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

PROTEIN_ACC_RE = re.compile(r"(?:XP|NP|YP|WP)_[0-9]+\.[0-9]+")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-fasta", required=True, type=Path, help="Input FASTA with original seq IDs.")
    p.add_argument(
        "--output-label-tsv",
        required=True,
        type=Path,
        help="Output TSV mapping seq_id to display_label.",
    )
    p.add_argument(
        "--output-itol-labels",
        required=True,
        type=Path,
        help="Output iTOL LABELS file.",
    )
    return p


def parse_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    continue
                ids.append(header.split(None, 1)[0])
    if not ids:
        raise ValueError(f"No FASTA headers found: {path}")
    return ids


def to_species_label(seq_id: str) -> str:
    protein_acc = ""
    left = seq_id
    m = PROTEIN_ACC_RE.search(seq_id)
    if m:
        protein_acc = m.group(0)
        left = seq_id[: m.start()]
    left = left.strip("_")
    tokens = [t for t in left.split("_") if t]

    if tokens and tokens[0] == "Poales":
        tokens = tokens[1:]

    if len(tokens) >= 2:
        species = f"{tokens[0]} {tokens[1]}"
        return f"{species} {protein_acc}".strip()

    cleaned = re.sub(r"_(?:CNL|CNLO|NL|NLR)$", "", seq_id)
    cleaned_label = cleaned.replace("_", " ").strip()
    return f"{cleaned_label} {protein_acc}".strip()


def write_label_tsv(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["seq_id", "display_label"])
        for seq_id, label in rows:
            w.writerow([seq_id, label])


def write_itol_labels(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fout:
        fout.write("LABELS\n")
        fout.write("SEPARATOR TAB\n")
        fout.write("DATA\n")
        for seq_id, label in rows:
            fout.write(f"{seq_id}\t{label}\n")


def main() -> None:
    args = build_parser().parse_args()
    seq_ids = parse_fasta_ids(args.input_fasta)
    rows = [(sid, to_species_label(sid)) for sid in seq_ids]

    write_label_tsv(args.output_label_tsv, rows)
    write_itol_labels(args.output_itol_labels, rows)

    print(f"Input IDs: {len(seq_ids)}")
    print(f"Label TSV: {args.output_label_tsv}")
    print(f"iTOL LABELS: {args.output_itol_labels}")


if __name__ == "__main__":
    main()
