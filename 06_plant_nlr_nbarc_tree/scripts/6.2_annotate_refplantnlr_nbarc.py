#!/usr/bin/env python3
"""Append NLRtracker domain architecture labels to RefPlantNLR NB-ARC IDs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path
from typing import Iterator


NBARC_COPY_PATTERN = re.compile(r"^(?P<base>.+?)(?:_\((?P<copy>[0-9]+)\))?$")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    seq_id: str | None = None
    sequence_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(sequence_parts)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header in {path}")
                seq_id = header.split(None, 1)[0]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
    if seq_id is not None:
        yield seq_id, "".join(sequence_parts)


def load_metadata(path: Path) -> dict[str, tuple[str, str]]:
    class_column = "Domain architecture simplified"
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"seqname", "Status", class_column}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {sorted(required)}")
        metadata: dict[str, tuple[str, str]] = {}
        for row in reader:
            seq_id = row["seqname"].strip()
            value = (row["Status"].strip(), row[class_column].strip())
            if seq_id in metadata and metadata[seq_id] != value:
                raise ValueError(f"Conflicting NLRtracker metadata for {seq_id}")
            metadata[seq_id] = value
    return metadata


def parse_nbarc_id(raw_id: str) -> tuple[str, int | None]:
    match = NBARC_COPY_PATTERN.fullmatch(raw_id)
    if not match:
        raise ValueError(f"Could not parse NB-ARC ID: {raw_id}")
    copy_text = match.group("copy")
    return match.group("base"), int(copy_text) if copy_text else None


def safe_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "", value)
    return label or "UNCLASSIFIED"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-nbarc-fasta", required=True, type=Path)
    parser.add_argument("--input-nlrtracker-tsv", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    metadata = load_metadata(args.input_nlrtracker_tsv)
    source_records = list(read_fasta(args.input_nbarc_fasta))
    parsed_records = [
        (raw_id, *parse_nbarc_id(raw_id), sequence)
        for raw_id, sequence in source_records
    ]
    base_counts = Counter(base for _, base, _, _ in parsed_records)
    encounter_counts: Counter[str] = Counter()
    output_ids: set[str] = set()
    retained = 0
    excluded_non_nlr = 0

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("w", encoding="utf-8") as output:
        for raw_id, base_id, source_copy, sequence in parsed_records:
            if base_id not in metadata:
                raise ValueError(f"No NLRtracker metadata for RefPlantNLR ID: {raw_id}")
            status, architecture_raw = metadata[base_id]
            if status != "NLR":
                excluded_non_nlr += 1
                continue

            encounter_counts[base_id] += 1
            architecture = safe_label(architecture_raw)
            if base_counts[base_id] > 1:
                copy_index = source_copy or encounter_counts[base_id]
                output_id = f"{base_id}_NBARC{copy_index}_{architecture}"
            else:
                output_id = f"{base_id}_{architecture}"

            if output_id in output_ids:
                raise ValueError(f"Duplicate annotated RefPlantNLR ID: {output_id}")
            output_ids.add(output_id)
            output.write(f">{output_id}\n{sequence}\n")
            retained += 1

    print(f"Source RefPlantNLR NB-ARC records: {len(source_records)}")
    print(f"Excluded records with Status other than NLR: {excluded_non_nlr}")
    print(f"Annotated RefPlantNLR NB-ARC records: {retained}")
    print(f"Output FASTA: {args.output_fasta}")


if __name__ == "__main__":
    main()
