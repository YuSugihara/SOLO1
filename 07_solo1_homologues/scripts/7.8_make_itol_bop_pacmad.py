#!/usr/bin/env python3
"""Build BOP/PACMAD iTOL annotation from sequence IDs using taxonkit."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tarfile
import urllib.request
from pathlib import Path

PROTEIN_ACC_RE = re.compile(r"(?:XP|NP|YP|WP)_[0-9]+\.[0-9]+")

DETAIL_CLADES_ORDER = [
    "PACMAD",
    "BOP:Triticodae",
    "BOP:Poodae",
    "BOP:other",
    "Poaceae_unresolved",
    "Non-Poaceae",
    "Unresolved",
]
DETAIL_CLADE_COLORS = {
    # High-contrast, colorblind-friendly palette for the 4 main classes.
    "PACMAD": "#CC79A7",
    "BOP:Triticodae": "#0072B2",
    "BOP:Poodae": "#E69F00",
    "BOP:other": "#009E73",
    "Poaceae_unresolved": "#FF7F00",
    "Non-Poaceae": "#33A02C",
    "Unresolved": "#BDBDBD",
}

BOP_MARKERS = {
    "BOP clade",
    "BEP clade",
    "Pooideae",
    "Oryzoideae",
    "Bambusoideae",
}
PACMAD_MARKERS = {
    "PACMAD clade",
    "Panicoideae",
    "Chloridoideae",
    "Arundinoideae",
    "Danthonioideae",
    "Aristidoideae",
    "Micrairoideae",
}
# Taxonomy naming can vary by NCBI updates; include common related ranks/tokens.
TRITICODAE_MARKERS = {"Triticodae", "Triticeae", "Triticinae"}
POODAE_MARKERS = {"Poodae", "Poeae", "Aveneae", "Bromeae", "Stipeae", "Meliceae"}
REQUIRED_DB_FILES = ("names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-fasta", required=True, type=Path, help="Input FASTA with sequence IDs.")
    p.add_argument(
        "--output-species-taxonomy-tsv",
        required=True,
        type=Path,
        help="Output TSV with species -> taxid/lineage/clade.",
    )
    p.add_argument(
        "--output-seqid-clade-tsv",
        required=True,
        type=Path,
        help="Output TSV with seq_id -> species -> clade mapping.",
    )
    p.add_argument(
        "--output-missing-tsv",
        required=True,
        type=Path,
        help="Output TSV with unresolved species/name lookup issues.",
    )
    p.add_argument(
        "--output-itol-colorstrip",
        required=True,
        type=Path,
        help="Output iTOL DATASET_COLORSTRIP file (with COLOR_BRANCHES=1).",
    )
    p.add_argument(
        "--dataset-label",
        required=False,
        default="BOP_PACMAD_clade",
        help="Dataset label in iTOL file.",
    )
    p.add_argument(
        "--taxonkit-data-dir",
        required=False,
        type=Path,
        help="Optional taxonkit --data-dir path.",
    )
    p.add_argument(
        "--auto-download-taxdump",
        action="store_true",
        help="If set, download and install NCBI taxdump when required files are missing.",
    )
    p.add_argument(
        "--taxdump-url",
        required=False,
        default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        help="URL to NCBI taxdump.tar.gz",
    )
    return p


def parse_fasta_ids(path: Path) -> list[str]:
    seq_ids: list[str] = []
    with path.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                if header:
                    seq_ids.append(header.split(None, 1)[0])
    if not seq_ids:
        raise ValueError(f"No FASTA headers found: {path}")
    return seq_ids


def seqid_to_species(seq_id: str) -> str:
    left = seq_id
    m = PROTEIN_ACC_RE.search(seq_id)
    if m:
        left = seq_id[: m.start()]
    left = left.strip("_")
    tokens = [t for t in left.split("_") if t]
    if tokens and tokens[0] == "Poales":
        tokens = tokens[1:]
    if len(tokens) >= 2:
        return f"{tokens[0]} {tokens[1]}"
    return ""


def run_taxonkit_name2taxid(species_names: list[str], data_dir: Path | None) -> dict[str, str]:
    cmd = ["taxonkit", "name2taxid"]
    if data_dir is not None:
        cmd.extend(["--data-dir", str(data_dir)])
    proc = subprocess.run(
        cmd,
        input="\n".join(species_names) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"taxonkit name2taxid failed:\n{proc.stderr.strip()}")

    out: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        name, taxid = parts[0].strip(), parts[1].strip()
        if name and name not in out:
            out[name] = taxid
    return out


def run_taxonkit_lineage(taxids: list[str], data_dir: Path | None) -> dict[str, str]:
    cmd = ["taxonkit", "lineage"]
    if data_dir is not None:
        cmd.extend(["--data-dir", str(data_dir)])
    proc = subprocess.run(
        cmd,
        input="\n".join(taxids) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"taxonkit lineage failed:\n{proc.stderr.strip()}")

    out: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        taxid, lineage = parts[0].strip(), parts[1].strip()
        if taxid and taxid not in out:
            out[taxid] = lineage
    return out


def has_taxonkit_db_files(data_dir: Path) -> bool:
    return all((data_dir / f).exists() for f in REQUIRED_DB_FILES)


def ensure_taxonkit_db(data_dir: Path, auto_download: bool, taxdump_url: str) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    if has_taxonkit_db_files(data_dir):
        return
    if not auto_download:
        missing = [f for f in REQUIRED_DB_FILES if not (data_dir / f).exists()]
        raise RuntimeError(
            "Taxonkit database files are missing in "
            f"{data_dir}. Missing: {', '.join(missing)}. "
            "Set --auto-download-taxdump to download them automatically."
        )

    tar_path = data_dir / "taxdump.tar.gz"
    urllib.request.urlretrieve(taxdump_url, tar_path)  # noqa: S310

    with tarfile.open(tar_path, "r:gz") as tf:
        members = {m.name: m for m in tf.getmembers()}
        for fname in REQUIRED_DB_FILES:
            if fname not in members:
                raise RuntimeError(f"Required file {fname} not found in downloaded taxdump.")
            src = tf.extractfile(members[fname])
            if src is None:
                raise RuntimeError(f"Failed to extract {fname} from taxdump.")
            with (data_dir / fname).open("wb") as fout:
                fout.write(src.read())

    if not has_taxonkit_db_files(data_dir):
        raise RuntimeError("Failed to prepare taxonkit database files after download.")


def tokenize_lineage(lineage: str) -> set[str]:
    return {tok.strip() for tok in lineage.split(";") if tok.strip()}


def classify_clade(lineage: str) -> tuple[str, str]:
    if not lineage:
        return ("Unresolved", "Unresolved")
    tokens = tokenize_lineage(lineage)
    if tokens & PACMAD_MARKERS:
        return ("PACMAD", "PACMAD")
    if tokens & BOP_MARKERS:
        if tokens & TRITICODAE_MARKERS:
            return ("BOP", "BOP:Triticodae")
        if tokens & POODAE_MARKERS:
            return ("BOP", "BOP:Poodae")
        return ("BOP", "BOP:other")
    if "Poaceae" in tokens:
        return ("Poaceae_unresolved", "Poaceae_unresolved")
    return ("Non-Poaceae", "Non-Poaceae")


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def write_itol_colorstrip(path: Path, dataset_label: str, seq_rows: list[list[str]]) -> None:
    present = {row[5] for row in seq_rows}
    legend_labels = [c for c in DETAIL_CLADES_ORDER if c in present] + sorted(
        c for c in present if c not in set(DETAIL_CLADES_ORDER)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fout:
        fout.write("DATASET_COLORSTRIP\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\t{dataset_label}\n")
        fout.write("COLOR\t#777777\n")
        fout.write("COLOR_BRANCHES\t1\n")
        fout.write("LEGEND_TITLE\tPoaceae superclades\n")
        fout.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_labels)) + "\n")
        fout.write("LEGEND_COLORS\t" + "\t".join(DETAIL_CLADE_COLORS[c] for c in legend_labels) + "\n")
        fout.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n")
        fout.write("DATA\n")
        for row in seq_rows:
            seq_id = row[0]
            detail_clade = row[5]
            color = row[6]
            fout.write(f"{seq_id}\t{color}\t{detail_clade}\n")


def main() -> None:
    args = build_parser().parse_args()
    data_dir = args.taxonkit_data_dir
    if args.auto_download_taxdump and data_dir is None:
        data_dir = Path.home() / ".taxonkit"

    if data_dir is not None:
        ensure_taxonkit_db(data_dir, args.auto_download_taxdump, args.taxdump_url)

    seq_ids = parse_fasta_ids(args.input_fasta)

    seq_to_species: dict[str, str] = {}
    for seq_id in seq_ids:
        seq_to_species[seq_id] = seqid_to_species(seq_id)

    unique_species = sorted({s for s in seq_to_species.values() if s})
    name_to_taxid = run_taxonkit_name2taxid(unique_species, data_dir)

    taxids = sorted({taxid for taxid in name_to_taxid.values() if taxid})
    taxid_to_lineage = run_taxonkit_lineage(taxids, data_dir)

    species_rows: list[list[str]] = []
    missing_rows: list[list[str]] = []
    species_to_info: dict[str, tuple[str, str, str, str]] = {}
    for species in unique_species:
        taxid = name_to_taxid.get(species, "")
        lineage = taxid_to_lineage.get(taxid, "") if taxid else ""
        clade_coarse, clade_detail = classify_clade(lineage)
        species_rows.append(
            [species, taxid, lineage, clade_coarse, clade_detail, DETAIL_CLADE_COLORS[clade_detail]]
        )
        species_to_info[species] = (taxid, lineage, clade_coarse, clade_detail)
        if not taxid:
            missing_rows.append([species, "name2taxid_not_found"])
        elif not lineage:
            missing_rows.append([species, "lineage_not_found"])

    seq_rows: list[list[str]] = []
    for seq_id in seq_ids:
        species = seq_to_species.get(seq_id, "")
        if not species:
            clade_coarse = "Unresolved"
            clade_detail = "Unresolved"
            taxid = ""
            lineage = ""
            note = "species_parse_failed"
            missing_rows.append([seq_id, note])
        else:
            taxid, lineage, clade_coarse, clade_detail = species_to_info.get(
                species, ("", "", "Unresolved", "Unresolved")
            )
            note = ""
        seq_rows.append(
            [
                seq_id,
                species,
                taxid,
                lineage,
                clade_coarse,
                clade_detail,
                DETAIL_CLADE_COLORS[clade_detail],
                note,
            ]
        )

    write_tsv(
        args.output_species_taxonomy_tsv,
        ["species", "taxid", "lineage", "clade_coarse", "clade_detail", "color_hex"],
        species_rows,
    )
    write_tsv(
        args.output_seqid_clade_tsv,
        ["seq_id", "species", "taxid", "lineage", "clade_coarse", "clade_detail", "color_hex", "note"],
        seq_rows,
    )
    write_tsv(args.output_missing_tsv, ["id_or_species", "issue"], sorted(missing_rows))
    write_itol_colorstrip(args.output_itol_colorstrip, args.dataset_label, seq_rows)

    counts_detail: dict[str, int] = {c: 0 for c in DETAIL_CLADES_ORDER}
    counts_coarse: dict[str, int] = {}
    for row in seq_rows:
        counts_detail[row[5]] = counts_detail.get(row[5], 0) + 1
        counts_coarse[row[4]] = counts_coarse.get(row[4], 0) + 1

    print(f"Input sequences: {len(seq_ids)}")
    print(f"Unique species: {len(unique_species)}")
    print(f"Species taxonomy TSV: {args.output_species_taxonomy_tsv}")
    print(f"SeqID clade TSV: {args.output_seqid_clade_tsv}")
    print(f"Missing TSV: {args.output_missing_tsv}")
    print(f"iTOL colorstrip (COLOR_BRANCHES=1): {args.output_itol_colorstrip}")
    if data_dir is not None:
        print(f"Taxonkit data-dir: {data_dir}")
        print("Taxonkit DB files checked/prepared: names.dmp,nodes.dmp,delnodes.dmp,merged.dmp")
    print("Coarse clade counts:")
    for c in sorted(counts_coarse):
        print(f"  {c}: {counts_coarse[c]}")
    print("Detailed clade counts:")
    for c in DETAIL_CLADES_ORDER:
        if counts_detail.get(c, 0) > 0:
            print(f"  {c}: {counts_detail[c]}")


if __name__ == "__main__":
    main()
