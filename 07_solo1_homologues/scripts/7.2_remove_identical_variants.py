#!/usr/bin/env python3
"""Remove redundant identical splicing variants from full-length NLR FASTA."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

XP_NUM_RE = re.compile(r"XP_(\d+)\.(\d+)")


@dataclass
class FastaRecord:
    index: int
    seq_id: str
    header: str
    seq: str


@dataclass
class Locus:
    assembly: str
    chrom: str
    strand: str
    start: int
    end: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-fasta", required=True, type=Path, help="Input full-length NLR FASTA")
    p.add_argument(
        "--seqid-assembly-tsv",
        required=True,
        type=Path,
        help="SEQ ID to assembly TSV (7.0_seqid_to_assembly.tsv; columns: seq_id, assembly_accession)",
    )
    p.add_argument(
        "--exon-boundaries-tsv",
        required=True,
        type=Path,
        help="Exon boundary TSV (7.0_exon_boundaries.tsv; columns include seq_id, seqid, start, end, strand, assembly_accession)",
    )
    p.add_argument("--output-fasta", required=True, type=Path, help="Output non-redundant FASTA")
    p.add_argument("--output-removed-tsv", required=True, type=Path, help="Removed redundant records")
    p.add_argument("--output-kept-tsv", required=True, type=Path, help="Kept records with locus info")
    return p.parse_args()


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header = ""
    seq_lines: list[str] = []
    idx = 0
    with path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                if header:
                    seq_id = header.split(None, 1)[0]
                    records.append(FastaRecord(index=idx, seq_id=seq_id, header=header, seq="".join(seq_lines)))
                    idx += 1
                header = line[1:].rstrip("\n\r")
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        seq_id = header.split(None, 1)[0]
        records.append(FastaRecord(index=idx, seq_id=seq_id, header=header, seq="".join(seq_lines)))
    return records


def load_seqid_to_assembly(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as fin:
        r = csv.DictReader(fin, delimiter="\t")
        required = {"seq_id", "assembly_accession"}
        if not r.fieldnames or not required.issubset(set(r.fieldnames)):
            raise ValueError("seqid-assembly TSV must include: seq_id, assembly_accession")
        for row in r:
            sid = (row.get("seq_id") or "").strip()
            asm = (row.get("assembly_accession") or "").strip()
            if sid:
                out[sid] = asm
    return out


def load_locus_info(path: Path) -> dict[str, Locus]:
    """
    Summarize exon coordinates to a single CDS span per seq_id.
    If multiple chromosome/strand groups are present, choose the group with largest covered span.
    """
    grouped: dict[str, dict[tuple[str, str, str], list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    with path.open("r", encoding="utf-8") as fin:
        r = csv.DictReader(fin, delimiter="\t")
        required = {"seq_id", "assembly_accession", "seqid", "start", "end", "strand"}
        if not r.fieldnames or not required.issubset(set(r.fieldnames)):
            raise ValueError(
                "exon-boundaries TSV must include: seq_id, assembly_accession, seqid, start, end, strand"
            )
        for row in r:
            sid = (row.get("seq_id") or "").strip()
            asm = (row.get("assembly_accession") or "").strip()
            chrom = (row.get("seqid") or "").strip()
            strand = (row.get("strand") or "").strip()
            if not sid or not asm or not chrom or not strand:
                continue
            try:
                s = int(row["start"])
                e = int(row["end"])
            except Exception:  # noqa: BLE001
                continue
            if s > e:
                s, e = e, s
            grouped[sid][(asm, chrom, strand)].append((s, e))

    locus_by_seq: dict[str, Locus] = {}
    for sid, locus_groups in grouped.items():
        best_key = None
        best_span = -1
        best_start = 0
        best_end = 0
        for key, intervals in locus_groups.items():
            starts = [x[0] for x in intervals]
            ends = [x[1] for x in intervals]
            span_start = min(starts)
            span_end = max(ends)
            span = span_end - span_start + 1
            if span > best_span:
                best_span = span
                best_key = key
                best_start = span_start
                best_end = span_end
        if best_key is None:
            continue
        asm, chrom, strand = best_key
        locus_by_seq[sid] = Locus(assembly=asm, chrom=chrom, strand=strand, start=best_start, end=best_end)
    return locus_by_seq


def overlaps(a: Locus, b: Locus) -> bool:
    if a.assembly != b.assembly:
        return False
    if a.chrom != b.chrom:
        return False
    if a.strand != b.strand:
        return False
    return not (a.end < b.start or b.end < a.start)


def connected_components(nodes: list[int], edges: dict[int, set[int]]) -> list[list[int]]:
    seen: set[int] = set()
    comps: list[list[int]] = []
    for n in nodes:
        if n in seen:
            continue
        stack = [n]
        seen.add(n)
        comp: list[int] = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in edges.get(cur, set()):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comps.append(sorted(comp))
    return comps


def representative_rank(record: FastaRecord) -> tuple[int, int, int]:
    """
    Smaller rank is preferred.
    Priority:
    1) smaller XP numeric part (XP_XXXXXXXX)
    2) smaller XP version
    3) earlier FASTA order
    """
    m = XP_NUM_RE.search(record.seq_id)
    if not m:
        return (10**18, 10**9, record.index)
    xp_num = int(m.group(1))
    xp_ver = int(m.group(2))
    return (xp_num, xp_ver, record.index)


def main() -> None:
    args = parse_args()
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_removed_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.output_kept_tsv.parent.mkdir(parents=True, exist_ok=True)

    records = read_fasta(args.input_fasta)
    seqid_to_assembly = load_seqid_to_assembly(args.seqid_assembly_tsv)
    locus_by_seq = load_locus_info(args.exon_boundaries_tsv)

    # Group by (assembly, amino-acid sequence). Only candidates with assembly info are considered.
    groups: dict[tuple[str, str], list[FastaRecord]] = defaultdict(list)
    keep_forced: set[str] = set()
    for r in records:
        asm = seqid_to_assembly.get(r.seq_id, "")
        if not asm:
            keep_forced.add(r.seq_id)
            continue
        groups[(asm, r.seq)].append(r)

    kept_ids: set[str] = set()
    removed_rows: list[list[str]] = []

    for (asm, aa_seq), recs in groups.items():
        if len(recs) == 1:
            kept_ids.add(recs[0].seq_id)
            continue

        nodes = list(range(len(recs)))
        edges: dict[int, set[int]] = defaultdict(set)
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                sid_i = recs[i].seq_id
                sid_j = recs[j].seq_id
                loc_i = locus_by_seq.get(sid_i)
                loc_j = locus_by_seq.get(sid_j)
                if not loc_i or not loc_j:
                    continue
                if overlaps(loc_i, loc_j):
                    edges[i].add(j)
                    edges[j].add(i)

        comps = connected_components(nodes, edges)
        for comp in comps:
            comp_records = [recs[i] for i in comp]
            rep = min(comp_records, key=representative_rank)
            kept_ids.add(rep.seq_id)

            # If component has no overlap edges (singleton) nothing is removed.
            if len(comp_records) == 1:
                continue

            for rr in comp_records:
                if rr.seq_id == rep.seq_id:
                    continue
                loc_r = locus_by_seq.get(rr.seq_id)
                loc_rep = locus_by_seq.get(rep.seq_id)
                removed_rows.append(
                    [
                        rr.seq_id,
                        rep.seq_id,
                        asm,
                        str(len(aa_seq)),
                        loc_r.chrom if loc_r else "",
                        loc_r.strand if loc_r else "",
                        str(loc_r.start) if loc_r else "",
                        str(loc_r.end) if loc_r else "",
                        loc_rep.chrom if loc_rep else "",
                        loc_rep.strand if loc_rep else "",
                        str(loc_rep.start) if loc_rep else "",
                        str(loc_rep.end) if loc_rep else "",
                        "same_assembly+identical_aa+genome_overlap",
                    ]
                )

    # Keep records without assembly mapping as-is.
    kept_ids.update(keep_forced)

    kept_records = [r for r in records if r.seq_id in kept_ids]
    with args.output_fasta.open("w", encoding="utf-8") as fout:
        for r in kept_records:
            fout.write(f">{r.header}\n")
            fout.write(f"{r.seq}\n")

    with args.output_removed_tsv.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(
            [
                "removed_seq_id",
                "representative_seq_id",
                "assembly_accession",
                "aa_length",
                "removed_chrom",
                "removed_strand",
                "removed_span_start",
                "removed_span_end",
                "rep_chrom",
                "rep_strand",
                "rep_span_start",
                "rep_span_end",
                "criterion",
            ]
        )
        w.writerows(sorted(removed_rows))

    with args.output_kept_tsv.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["seq_id", "assembly_accession", "chrom", "strand", "span_start", "span_end", "aa_length"])
        for r in kept_records:
            asm = seqid_to_assembly.get(r.seq_id, "")
            loc = locus_by_seq.get(r.seq_id)
            w.writerow(
                [
                    r.seq_id,
                    asm,
                    loc.chrom if loc else "",
                    loc.strand if loc else "",
                    str(loc.start) if loc else "",
                    str(loc.end) if loc else "",
                    len(r.seq),
                ]
            )

    removed_count = len(records) - len(kept_records)
    print(f"Input records: {len(records)}")
    print(f"Kept records: {len(kept_records)}")
    print(f"Removed records: {removed_count}")
    print(f"Output FASTA: {args.output_fasta}")
    print(f"Removed TSV: {args.output_removed_tsv}")
    print(f"Kept TSV: {args.output_kept_tsv}")


if __name__ == "__main__":
    main()
