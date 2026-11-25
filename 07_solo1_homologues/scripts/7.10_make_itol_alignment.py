#!/usr/bin/env python3
"""Build an iTOL DATASET_ALIGNMENT file from an aligned FASTA."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input-alignment-fasta",
        required=True,
        type=Path,
        help="Aligned FASTA to visualize in iTOL (tip IDs must match tree leaf IDs).",
    )
    p.add_argument(
        "--output-itol-alignment",
        required=True,
        type=Path,
        help="Output iTOL DATASET_ALIGNMENT file.",
    )
    p.add_argument(
        "--output-summary-tsv",
        required=False,
        type=Path,
        help="Optional summary TSV (record count, alignment length stats).",
    )
    p.add_argument(
        "--dataset-label",
        required=False,
        default="full_length_alignment",
        help="Dataset label shown in iTOL.",
    )
    p.add_argument(
        "--color",
        required=False,
        default="#3B6FB6",
        help="Dataset color shown in iTOL.",
    )
    p.add_argument(
        "--consensus-threshold",
        required=False,
        type=int,
        default=50,
        help="Consensus threshold (0-100).",
    )
    p.add_argument(
        "--ignore-gaps",
        required=False,
        type=int,
        choices=[0, 1],
        default=1,
        help="Set to 1 to ignore gaps in consensus calculation.",
    )
    p.add_argument(
        "--allow-variable-length",
        action="store_true",
        help="Allow variable sequence lengths. Default: require equal lengths.",
    )
    return p


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    seq_lines: list[str] = []
    with path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                if header:
                    seq_id = header.split(None, 1)[0]
                    records.append((seq_id, "".join(seq_lines)))
                header = line[1:].rstrip("\n\r")
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        seq_id = header.split(None, 1)[0]
        records.append((seq_id, "".join(seq_lines)))
    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    return records


def write_itol_alignment(
    path: Path,
    records: list[tuple[str, str]],
    dataset_label: str,
    color: str,
    consensus_threshold: int,
    ignore_gaps: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fout:
        fout.write("DATASET_ALIGNMENT\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\t{dataset_label}\n")
        fout.write(f"COLOR\t{color}\n")
        fout.write(f"CONSENSUS_THRESHOLD\t{consensus_threshold}\n")
        fout.write(f"IGNORE_GAPS\t{ignore_gaps}\n")
        fout.write("DATA\n")
        for seq_id, seq in records:
            fout.write(f">{seq_id}\n{seq}\n")


def write_summary(path: Path, lengths: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    uniform = "yes" if len(set(lengths)) == 1 else "no"
    with path.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["record_count", len(lengths)])
        w.writerow(["alignment_length_min", min(lengths)])
        w.writerow(["alignment_length_max", max(lengths)])
        w.writerow(["alignment_length_uniform", uniform])


def main() -> None:
    args = build_parser().parse_args()
    if not (0 <= args.consensus_threshold <= 100):
        raise ValueError("--consensus-threshold must be between 0 and 100.")

    records = parse_fasta(args.input_alignment_fasta)
    lengths = [len(seq) for _, seq in records]
    uniform = len(set(lengths)) == 1
    if not uniform and not args.allow_variable_length:
        raise ValueError(
            "Alignment sequences are not all the same length. "
            "Use --allow-variable-length to bypass this check."
        )

    write_itol_alignment(
        args.output_itol_alignment,
        records,
        args.dataset_label,
        args.color,
        args.consensus_threshold,
        args.ignore_gaps,
    )
    if args.output_summary_tsv:
        write_summary(args.output_summary_tsv, lengths)

    print(f"Input records: {len(records)}")
    print(f"Alignment length min/max: {min(lengths)}/{max(lengths)}")
    print(f"Uniform alignment length: {'yes' if uniform else 'no'}")
    print(f"iTOL alignment dataset: {args.output_itol_alignment}")
    if args.output_summary_tsv:
        print(f"Summary TSV: {args.output_summary_tsv}")


if __name__ == "__main__":
    main()
