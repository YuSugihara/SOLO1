#!/usr/bin/env python3
"""Build exon/domain summary table and iTOL annotation inputs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

PROTEIN_ACC_RE = re.compile(r"(?:XP|NP|YP|WP)_[0-9]+\.[0-9]+")
DOMAIN_KEYS = ("CC", "NBARC", "LRR")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target-fasta", required=True, type=Path, help="Target FASTA (tree tips)")
    p.add_argument(
        "--seqid-assembly-tsv",
        required=True,
        type=Path,
        help="SEQ ID to assembly TSV (7.0_seqid_to_assembly.tsv)",
    )
    p.add_argument(
        "--exon-boundaries-tsv",
        required=True,
        type=Path,
        help="Exon boundaries TSV (7.0_exon_boundaries.tsv)",
    )
    p.add_argument(
        "--nlrtracker-root",
        required=True,
        type=Path,
        help=(
            "Directory holding the NLRtracker domain tables, one file per "
            "assembly, named nlrtracker_<accession>_protein.faa_Domains.tsv."
        ),
    )
    p.add_argument(
        "--alignment-fasta",
        required=True,
        type=Path,
        help="Aligned FASTA used for anchor-based LRR length measurements",
    )
    p.add_argument(
        "--reference-seq-id",
        required=False,
        default="Poales_Aegilops_tauschii_subsp._strangulata_XP_020170543.1_CNLO",
        help="Reference sequence ID used to define anchor coordinates on alignment",
    )
    p.add_argument(
        "--reference-lrr-start",
        required=False,
        type=int,
        default=510,
        help="1-based LRR start position on ungapped reference sequence",
    )
    p.add_argument(
        "--reference-lrr-end",
        required=False,
        type=int,
        default=833,
        help="1-based LRR end position on ungapped reference sequence",
    )
    p.add_argument("--output-table-tsv", required=True, type=Path, help="Output summary table TSV")
    p.add_argument("--output-missing-tsv", required=True, type=Path, help="Output missing/info TSV")
    p.add_argument(
        "--output-lrr-tail-fasta",
        required=False,
        type=Path,
        help="Optional output FASTA of ungapped sequences after reference LRR end (exclusive)",
    )
    p.add_argument(
        "--output-itol-prefix",
        required=False,
        type=Path,
        help="Optional output prefix for iTOL simplebar/multibar file",
    )
    return p.parse_args()


def read_fasta_lengths(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    header = ""
    seq_lines: list[str] = []
    with path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                if header:
                    sid = header.split(None, 1)[0]
                    out[sid] = len("".join(seq_lines))
                header = line[1:].rstrip("\n\r")
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        sid = header.split(None, 1)[0]
        out[sid] = len("".join(seq_lines))
    return out


def read_alignment(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    header = ""
    seq_lines: list[str] = []
    with path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                if header:
                    sid = header.split(None, 1)[0]
                    out[sid] = "".join(seq_lines)
                header = line[1:].rstrip("\n\r")
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        sid = header.split(None, 1)[0]
        out[sid] = "".join(seq_lines)
    return out


def ref_pos_to_aln_col(ref_aln_seq: str, ref_pos_1based: int) -> int:
    if ref_pos_1based <= 0:
        raise ValueError("Reference position must be >= 1")
    aa_count = 0
    for i, aa in enumerate(ref_aln_seq):
        if aa != "-":
            aa_count += 1
            if aa_count == ref_pos_1based:
                return i
    raise ValueError(f"Reference position {ref_pos_1based} exceeds ungapped reference length {aa_count}")


def extract_protein_acc(seq_id: str) -> str:
    m = PROTEIN_ACC_RE.search(seq_id)
    if not m:
        raise ValueError(f"No protein accession in SEQ ID: {seq_id}")
    return m.group(0)


def load_seqid_assembly(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as fin:
        r = csv.DictReader(fin, delimiter="\t")
        required = {"seq_id", "assembly_accession"}
        if not r.fieldnames or not required.issubset(set(r.fieldnames)):
            raise ValueError("seqid-assembly TSV must include seq_id and assembly_accession")
        for row in r:
            sid = (row.get("seq_id") or "").strip()
            asm = (row.get("assembly_accession") or "").strip()
            if sid:
                out[sid] = asm
    return out


def load_exon_counts(path: Path) -> dict[str, int]:
    by_seq: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as fin:
        r = csv.DictReader(fin, delimiter="\t")
        required = {"seq_id", "exon_number_transcript"}
        if not r.fieldnames or not required.issubset(set(r.fieldnames)):
            raise ValueError("exon-boundaries TSV must include seq_id and exon_number_transcript")
        for row in r:
            sid = (row.get("seq_id") or "").strip()
            ex = (row.get("exon_number_transcript") or "").strip()
            if sid and ex:
                by_seq[sid].add(ex)
    return {k: len(v) for k, v in by_seq.items()}


def find_domain_tsv(nlrtracker_root: Path, assembly: str) -> Path | None:
    direct = (
        nlrtracker_root
        / f"nlrtracker_{assembly}_protein.faa"
        / f"nlrtracker_{assembly}_protein.faa_Domains.tsv"
    )
    if direct.exists():
        return direct
    candidates = list(nlrtracker_root.glob(f"**/nlrtracker_{assembly}_protein.faa_Domains.tsv"))
    if candidates:
        return candidates[0]
    return None


def load_domains_for_assembly(domain_tsv: Path) -> dict[str, dict[str, list[tuple[int, int]]]]:
    # protein -> domain -> intervals
    out: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    with domain_tsv.open("r", encoding="utf-8") as fin:
        r = csv.DictReader(fin, delimiter="\t")
        required = {"seqname", "description", "start", "end"}
        if not r.fieldnames or not required.issubset(set(r.fieldnames)):
            raise ValueError(f"Domain TSV missing required columns: {domain_tsv}")
        for row in r:
            protein = (row.get("seqname") or "").strip()
            desc = (row.get("description") or "").strip().upper()
            if not protein or desc not in DOMAIN_KEYS:
                continue
            try:
                s = int(row["start"])
                e = int(row["end"])
            except Exception:  # noqa: BLE001
                continue
            if s > e:
                s, e = e, s
            out[protein][desc].append((s, e))
    return out


def merged_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[tuple[int, int]] = []
    for s, e in sorted(set(intervals)):
        if not merged or s > merged[-1][1] + 1:
            merged.append((s, e))
        else:
            ps, pe = merged[-1]
            merged[-1] = (ps, max(pe, e))
    return sum(e - s + 1 for s, e in merged)


def write_itol_simplebar(path: Path, rows: list[dict[str, str]], field: str, label: str, color: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fout:
        fout.write("DATASET_SIMPLEBAR\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\t{label}\n")
        fout.write(f"COLOR\t{color}\n")
        fout.write("DATA\n")
        for row in rows:
            fout.write(f"{row['seq_id']}\t{row[field]}\n")


def write_itol_multibar(
    path: Path,
    rows: list[dict[str, str]],
    fields: list[tuple[str, str, str]],
    dataset_label: str,
    dataset_color: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fout:
        fout.write("DATASET_MULTIBAR\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\t{dataset_label}\n")
        fout.write(f"COLOR\t{dataset_color}\n")
        fout.write("FIELD_LABELS\t" + "\t".join(lbl for _, lbl, _ in fields) + "\n")
        fout.write("FIELD_COLORS\t" + "\t".join(col for _, _, col in fields) + "\n")
        fout.write("DATA\n")
        for row in rows:
            vals = [row[field] for field, _, _ in fields]
            fout.write("\t".join([row["seq_id"], *vals]) + "\n")


def write_fasta_no_wrap(path: Path, seqs: dict[str, str], order: list[str]) -> tuple[int, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    skipped_empty = 0
    with path.open("w", encoding="utf-8") as fout:
        for sid in order:
            seq = seqs.get(sid, "")
            if not seq:
                skipped_empty += 1
                continue
            fout.write(f">{sid}\n{seq}\n")
            written += 1
    return written, skipped_empty


def main() -> None:
    args = parse_args()
    args.output_table_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.output_missing_tsv.parent.mkdir(parents=True, exist_ok=True)
    if args.output_itol_prefix is not None:
        args.output_itol_prefix.parent.mkdir(parents=True, exist_ok=True)

    fasta_lengths = read_fasta_lengths(args.target_fasta)
    alignment = read_alignment(args.alignment_fasta)
    seqid_to_assembly = load_seqid_assembly(args.seqid_assembly_tsv)
    exon_counts = load_exon_counts(args.exon_boundaries_tsv)

    if args.reference_seq_id not in alignment:
        raise ValueError(
            f"Reference seq ID not found in alignment: {args.reference_seq_id} ({args.alignment_fasta})"
        )
    ref_aln = alignment[args.reference_seq_id]
    ref_lrr_start_col = ref_pos_to_aln_col(ref_aln, args.reference_lrr_start)
    ref_lrr_end_col = ref_pos_to_aln_col(ref_aln, args.reference_lrr_end)
    if ref_lrr_start_col > ref_lrr_end_col:
        ref_lrr_start_col, ref_lrr_end_col = ref_lrr_end_col, ref_lrr_start_col

    domain_cache: dict[str, tuple[Path | None, dict[str, dict[str, list[tuple[int, int]]]]]] = {}
    lrr_tail_after_end: dict[str, str] = {}

    rows: list[dict[str, str]] = []
    missing_rows: list[list[str]] = []

    for seq_id, prot_len in fasta_lengths.items():
        try:
            protein = extract_protein_acc(seq_id)
        except Exception as exc:  # noqa: BLE001
            missing_rows.append([seq_id, "", "", f"protein_accession_parse_error:{exc}"])
            continue

        assembly = seqid_to_assembly.get(seq_id, "")
        if not assembly:
            missing_rows.append([seq_id, protein, "", "assembly_not_found"])
            continue

        if assembly not in domain_cache:
            dom_path = find_domain_tsv(args.nlrtracker_root, assembly)
            if dom_path is None:
                domain_cache[assembly] = (None, {})
            else:
                domain_cache[assembly] = (dom_path, load_domains_for_assembly(dom_path))

        dom_path, dom_map = domain_cache[assembly]
        protein_domains = dom_map.get(protein, {})

        cc_intervals = protein_domains.get("CC", [])
        nb_intervals = protein_domains.get("NBARC", [])
        lrr_intervals = protein_domains.get("LRR", [])

        cc_len = merged_length(cc_intervals)
        nb_len = merged_length(nb_intervals)
        lrr_len_domain = merged_length(lrr_intervals)

        lrr_start = min((s for s, _ in lrr_intervals), default=0)
        lrr_end = max((e for _, e in lrr_intervals), default=0)
        aln_seq = alignment.get(seq_id, "")
        if aln_seq:
            if len(aln_seq) != len(ref_aln):
                missing_rows.append([seq_id, protein, assembly, "alignment_length_mismatch_with_reference"])
                lrr_len_aln = 0
                lrr_end_to_end = 0
                lrr_tail_after_end[seq_id] = ""
            else:
                lrr_region = aln_seq[ref_lrr_start_col : ref_lrr_end_col + 1]
                lrr_len_aln = sum(1 for aa in lrr_region if aa != "-")
                tail_after_end_region = aln_seq[ref_lrr_end_col + 1 :]
                tail_seq = "".join(aa for aa in tail_after_end_region if aa != "-")
                lrr_end_to_end = len(tail_seq)
                lrr_tail_after_end[seq_id] = tail_seq
        else:
            lrr_len_aln = 0
            lrr_end_to_end = 0
            lrr_tail_after_end[seq_id] = ""
            missing_rows.append([seq_id, protein, assembly, "sequence_not_found_in_alignment"])

        exon_count = exon_counts.get(seq_id, 0)

        rows.append(
            {
                "seq_id": seq_id,
                "protein_accession": protein,
                "assembly_accession": assembly,
                "protein_length_aa": str(prot_len),
                "exon_count": str(exon_count),
                "cc_length_aa": str(cc_len),
                "nbarc_length_aa": str(nb_len),
                "lrr_length_aa": str(lrr_len_aln),
                "lrr_length_domain_aa": str(lrr_len_domain),
                "lrr_start": str(lrr_start) if lrr_start else "",
                "lrr_end": str(lrr_end) if lrr_end else "",
                "lrr_end_to_end_length_aa": str(lrr_end_to_end),
            }
        )

        if dom_path is None:
            missing_rows.append([seq_id, protein, assembly, "domain_tsv_not_found"])
        elif protein not in dom_map:
            missing_rows.append([seq_id, protein, assembly, "protein_not_found_in_domain_tsv"])
        elif not lrr_intervals:
            missing_rows.append([seq_id, protein, assembly, "lrr_not_found"])

    rows_sorted = sorted(rows, key=lambda x: x["seq_id"])

    with args.output_table_tsv.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(
            [
                "seq_id",
                "protein_accession",
                "assembly_accession",
                "protein_length_aa",
                "exon_count",
                "cc_length_aa",
                "nbarc_length_aa",
                "lrr_length_aa",
                "lrr_length_domain_aa",
                "lrr_start",
                "lrr_end",
                "lrr_end_to_end_length_aa",
                "lrr_anchor_reference_seq_id",
                "lrr_anchor_reference_start",
                "lrr_anchor_reference_end",
            ]
        )
        for r in rows_sorted:
            w.writerow(
                [
                    r["seq_id"],
                    r["protein_accession"],
                    r["assembly_accession"],
                    r["protein_length_aa"],
                    r["exon_count"],
                    r["cc_length_aa"],
                    r["nbarc_length_aa"],
                    r["lrr_length_aa"],
                    r["lrr_length_domain_aa"],
                    r["lrr_start"],
                    r["lrr_end"],
                    r["lrr_end_to_end_length_aa"],
                    args.reference_seq_id,
                    str(args.reference_lrr_start),
                    str(args.reference_lrr_end),
                ]
            )

    with args.output_missing_tsv.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["seq_id", "protein_accession", "assembly_accession", "note"])
        for row in sorted(missing_rows):
            w.writerow(row)

    generated = []

    if args.output_itol_prefix is not None:
        stacked_file = Path(f"{args.output_itol_prefix}.domain_lengths_stacked.txt")
        write_itol_multibar(
            stacked_file,
            rows_sorted,
            fields=[
                ("cc_length_aa", "CC_len", "#3B6FB6"),
                ("nbarc_length_aa", "NBARC_len", "#D68C45"),
                ("lrr_length_aa", "LRR_len", "#4A9A8F"),
                ("lrr_end_to_end_length_aa", "LRR_end_to_Cterm_len", "#7E6BA8"),
            ],
            dataset_label="domain_lengths_stacked",
            dataset_color="#5B7A9A",
        )
        generated.append(stacked_file.as_posix())

    if args.output_lrr_tail_fasta:
        order = [r["seq_id"] for r in rows_sorted]
        written, skipped_empty = write_fasta_no_wrap(args.output_lrr_tail_fasta, lrr_tail_after_end, order)
    else:
        written, skipped_empty = 0, 0

    print(f"Target sequences: {len(fasta_lengths)}")
    print(
        "LRR anchor on reference: "
        f"{args.reference_seq_id} positions {args.reference_lrr_start}-{args.reference_lrr_end}"
    )
    print(f"Summary table: {args.output_table_tsv}")
    print(f"Missing/info table: {args.output_missing_tsv}")
    if args.output_lrr_tail_fasta:
        print(f"LRR-tail FASTA (after end, exclusive): {args.output_lrr_tail_fasta}")
        print(f"LRR-tail FASTA records written: {written}")
        print(f"LRR-tail FASTA records skipped (empty tail): {skipped_empty}")
    if generated:
        print("iTOL separated files:")
        for p in generated:
            print(f"  - {p}")

    print(
        f"Step 7.7 build_domain_length_table completed: sequences={len(rows_sorted)}, "
        f"missing={len(missing_rows)}, ref_lrr_cols={ref_lrr_start_col + 1}-{ref_lrr_end_col + 1}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
