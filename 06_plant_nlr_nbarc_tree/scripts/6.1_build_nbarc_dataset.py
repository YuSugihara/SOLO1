#!/usr/bin/env python3
"""Step 6.1 script: build full-length NLR and NB-ARC FASTAs for the selected genomes.

For each assembly listed in the manifest, retain the proteins that NLRtracker
classifies with `Status = NLR` and write their full-length and NB-ARC
sequences. Identifiers carry the taxon prefix and end in the simplified
domain architecture, so that the architecture can be read back from a tree
tip; a protein with several NB-ARC regions gets one record per region.
"""

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


def load_selected_manifest(manifest: Path) -> list[dict[str, str]]:
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "assembly_accession",
            "species",
            "fasta_taxon_prefix",
            "status_nlr_count",
            "role",
            "include",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Manifest is missing required columns: "
                f"{sorted(required - set(reader.fieldnames or []))}"
            )
        rows = [
            row
            for row in reader
            if row["include"].strip().lower() in {"yes", "true", "1"}
        ]
    if not rows:
        raise ValueError("No manifest rows have include=yes.")
    return rows


def load_tracker_metadata(path: Path) -> dict[str, tuple[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        class_column = "Domain architecture simplified"
        required = {"seqname", "Status", class_column}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must contain columns: {sorted(required)}"
            )
        metadata: dict[str, tuple[str, str]] = {}
        for row in reader:
            seq_id = row["seqname"].strip()
            value = (row["Status"].strip(), row[class_column].strip())
            if seq_id in metadata and metadata[seq_id] != value:
                raise ValueError(f"Conflicting NLRtracker metadata for {seq_id}")
            metadata[seq_id] = value
    return metadata


def safe_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "", value)
    return label or "UNCLASSIFIED"


def parse_nbarc_id(raw_id: str) -> tuple[str, int | None]:
    match = NBARC_COPY_PATTERN.fullmatch(raw_id)
    if not match:
        raise ValueError(f"Could not parse NB-ARC ID: {raw_id}")
    copy_text = match.group("copy")
    return match.group("base"), int(copy_text) if copy_text else None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-manifest", required=True, type=Path)
    parser.add_argument(
        "--nlrtracker-dir",
        required=True,
        type=Path,
        help=(
            "Directory holding one subdirectory per assembly accession, each "
            "with NLR.fasta, NBARC.fasta and NLRtracker.tsv."
        ),
    )
    parser.add_argument("--output-full-length-fasta", required=True, type=Path)
    parser.add_argument("--output-nbarc-fasta", required=True, type=Path)
    parser.add_argument("--output-summary-tsv", required=True, type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rows = load_selected_manifest(args.genome_manifest)
    args.output_full_length_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_nbarc_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str | int]] = []
    output_full_ids: set[str] = set()
    output_nbarc_ids: set[str] = set()
    total_full = 0
    total_nbarc = 0
    total_excluded_non_nlr = 0

    with args.output_full_length_fasta.open(
        "w", encoding="utf-8"
    ) as full_output, args.output_nbarc_fasta.open(
        "w", encoding="utf-8"
    ) as nbarc_output:
        for row in rows:
            assembly = row["assembly_accession"].strip()
            prefix = row["fasta_taxon_prefix"].strip()
            input_dir = args.nlrtracker_dir / assembly
            full_path = input_dir / "NLR.fasta"
            nbarc_path = input_dir / "NBARC.fasta"
            tracker_path = input_dir / "NLRtracker.tsv"

            for required_path in (full_path, nbarc_path, tracker_path):
                if not required_path.is_file() or required_path.stat().st_size == 0:
                    raise FileNotFoundError(
                        f"Required staged input is missing or empty: {required_path}"
                    )

            metadata = load_tracker_metadata(tracker_path)
            source_full_records = list(read_fasta(full_path))
            source_full_accessions = {
                seq_id for seq_id, _ in source_full_records
            }
            if len(source_full_accessions) != len(source_full_records):
                raise ValueError(f"Duplicate full-length IDs in {full_path}")

            for accession in source_full_accessions:
                if accession not in metadata:
                    raise ValueError(
                        f"{assembly}: no NLRtracker metadata for {accession}"
                    )

            full_records = [
                (accession, sequence)
                for accession, sequence in source_full_records
                if metadata[accession][0] == "NLR"
            ]
            full_accessions = {seq_id for seq_id, _ in full_records}
            excluded_non_nlr = len(source_full_records) - len(full_records)

            expected_full = int(row["status_nlr_count"])
            if len(full_records) != expected_full:
                raise ValueError(
                    f"{assembly}: expected {expected_full} Status=NLR records "
                    f"from the manifest, observed {len(full_records)}"
                )

            for accession, sequence in full_records:
                output_id = (
                    f"{prefix}_{accession}_{safe_label(metadata[accession][1])}"
                )
                if output_id in output_full_ids:
                    raise ValueError(f"Duplicate combined full-length ID: {output_id}")
                output_full_ids.add(output_id)
                full_output.write(f">{output_id}\n{sequence}\n")

            nbarc_records = list(read_fasta(nbarc_path))
            parsed_nbarc = [
                (raw_id, *parse_nbarc_id(raw_id), sequence)
                for raw_id, sequence in nbarc_records
            ]
            selected_nbarc = [
                record for record in parsed_nbarc if record[1] in full_accessions
            ]
            copy_counts = Counter(base for _, base, _, _ in selected_nbarc)
            encounter_counts: Counter[str] = Counter()
            observed_nbarc_accessions: set[str] = set()

            for raw_id, accession, _, _ in parsed_nbarc:
                if accession not in source_full_accessions:
                    raise ValueError(
                        f"{assembly}: NB-ARC {raw_id} has no matching full-length NLR"
                    )

            for raw_id, accession, source_copy, sequence in selected_nbarc:
                observed_nbarc_accessions.add(accession)
                encounter_counts[accession] += 1
                architecture = safe_label(metadata[accession][1])
                if copy_counts[accession] > 1:
                    copy_index = source_copy or encounter_counts[accession]
                    output_id = (
                        f"{prefix}_{accession}_NBARC{copy_index}_{architecture}"
                    )
                else:
                    output_id = f"{prefix}_{accession}_{architecture}"
                if output_id in output_nbarc_ids:
                    raise ValueError(f"Duplicate combined NB-ARC ID: {output_id}")
                output_nbarc_ids.add(output_id)
                nbarc_output.write(f">{output_id}\n{sequence}\n")

            without_nbarc = len(full_accessions - observed_nbarc_accessions)
            summary_rows.append(
                {
                    "assembly_accession": assembly,
                    "species": row["species"],
                    "role": row["role"],
                    "source_nlr_and_degenerate_count": len(source_full_records),
                    "excluded_non_nlr_count": excluded_non_nlr,
                    "retained_status_nlr_count": len(full_records),
                    "retained_nbarc_record_count": len(selected_nbarc),
                    "full_length_ids_without_nbarc": without_nbarc,
                }
            )
            total_full += len(full_records)
            total_nbarc += len(selected_nbarc)
            total_excluded_non_nlr += excluded_non_nlr

    fields = [
        "assembly_accession",
        "species",
        "role",
        "source_nlr_and_degenerate_count",
        "excluded_non_nlr_count",
        "retained_status_nlr_count",
        "retained_nbarc_record_count",
        "full_length_ids_without_nbarc",
    ]
    with args.output_summary_tsv.open(
        "w", encoding="utf-8", newline=""
    ) as summary_output:
        writer = csv.DictWriter(
            summary_output,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Selected assemblies: {len(rows)}")
    print(f"Excluded non-NLR full-length records: {total_excluded_non_nlr}")
    print(f"Combined full-length NLR records: {total_full}")
    print(f"Combined NB-ARC records: {total_nbarc}")
    print(f"Full-length FASTA: {args.output_full_length_fasta}")
    print(f"NB-ARC FASTA: {args.output_nbarc_fasta}")
    print(f"Summary TSV: {args.output_summary_tsv}")


if __name__ == "__main__":
    main()
