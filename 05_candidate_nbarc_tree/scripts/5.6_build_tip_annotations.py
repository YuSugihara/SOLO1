#!/usr/bin/env python3
"""Step 5.6 script: compile the tip annotation table for the candidate tree.

Builds a unified tab-separated annotation table keyed on the 296 tip labels of
the midpoint-rooted NB-ARC phylogeny, combining structural ranking metrics,
sensor-clade assignments, SOLO1 clade membership, taxonomic plant order, and
functional cell-death screening outcomes.
"""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

TOP_ORDERS = ("Poales", "Solanales", "Lamiales", "Caryophyllales")
ORDER_OTHER_LABEL = "Other"
CELL_DEATH_STAGE = "structural ranking"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--tree",
        required=True,
        type=Path,
        help="Midpoint-rooted tree whose tip labels the annotations key on.",
    )
    p.add_argument(
        "--table",
        required=True,
        type=Path,
        help="Ranked and annotated table from Step 4.1.",
    )
    p.add_argument(
        "--plant-order",
        required=True,
        type=Path,
        help="Table assigning each species to its taxonomic order.",
    )
    p.add_argument(
        "--solo1-ids",
        required=True,
        type=Path,
        help="List of SOLO1 clade accessions (Step 7.0).",
    )
    p.add_argument(
        "--output-tsv",
        required=True,
        type=Path,
        help="Path to write the unified tip annotation table (.tsv).",
    )
    return p.parse_args()


def fail(message):
    sys.exit(f"error: {message}")


def read_tip_labels(path):
    """Return the tip labels of a Newick tree, in file order."""
    text = path.read_text(encoding="utf-8")
    labels, token, in_label = [], [], False
    for char in text:
        if char in "(,":
            token, in_label = [], True
        elif char in ":)":
            if in_label and token:
                labels.append("".join(token).strip())
            token, in_label = [], False
        elif in_label:
            token.append(char)
    if not labels:
        fail(f"no tip labels found in {path}")
    if len(set(labels)) != len(labels):
        fail(f"{path} has duplicate tip labels")
    return labels


def read_table(path):
    rows = {}
    needed = (
        "protein_ID",
        "weighted_rank",
        "sensor_clade",
        "screening_class",
        "D_APEX",
        "ipTM_LCB",
        "Cell_death",
        "WT_domesticated_autoactivity",
        "screening_stage",
        "is_reference",
    )
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for column in needed:
            if column not in header:
                fail(f"{path} has no '{column}' column")
        index = {c: header.index(c) for c in needed}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            row = {c: fields[i] for c, i in index.items()}
            if row["is_reference"] == "True":
                continue
            rows[row["protein_ID"]] = row
    return rows


def read_plant_order(path):
    """Return the taxonomic order of each species, keyed on the species name."""
    orders = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for column in ("species", "order"):
            if column not in header:
                fail(f"{path} has no '{column}' column")
        species_at, order_at = header.index("species"), header.index("order")
        for number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(species_at, order_at):
                fail(f"{path} line {number} has too few fields")
            species = fields[species_at].strip()
            if species in orders:
                fail(f"{path} line {number} repeats the species {species}")
            orders[species] = fields[order_at].strip()
    if not orders:
        fail(f"{path} holds no species")
    return orders


def read_solo1_ids(path):
    """Extract SOLO1 clade protein accessions from list file."""
    accs = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            m = re.search(r"[XN]P_\d+\.\d+", line)
            if m:
                accs.add(m.group(0))
    if not accs:
        fail(f"no accessions found in {path}")
    return accs


def main():
    args = parse_args()
    labels = read_tip_labels(args.tree)
    table = read_table(args.table)
    plant_order = read_plant_order(args.plant_order)
    solo1_accs = read_solo1_ids(args.solo1_ids)

    columns = (
        "tip_label",
        "species",
        "assembly",
        "protein_ID",
        "plant_order",
        "plant_order_group",
        "sensor_clade",
        "is_solo1_clade",
        "weighted_rank",
        "ipTM_LCB",
        "D_APEX",
        "screening_stage",
        "Cell_death",
        "WT_domesticated_autoactivity",
        "screening_class",
    )

    records = []
    order_group_counts = Counter()
    all_orders = set()
    solo1_tip_count = 0
    tested_counts = Counter()

    for label in labels:
        parts = label.split("|")
        if len(parts) != 3:
            fail(f"tip label {label!r} is not <species>|<assembly>|<protein>")
        species, assembly, protein = parts
        if protein not in table:
            fail(f"{protein} is a tip of the tree but is absent from {args.table}")
        if species not in plant_order:
            fail(
                f"{species} is a tip of the tree but is absent from "
                f"{args.plant_order}"
            )
        row = table[protein]
        if row["screening_stage"] == CELL_DEATH_STAGE:
            if row["Cell_death"] not in ("Yes", "No"):
                fail(
                    f"{label} is on the {CELL_DEATH_STAGE} stage "
                    f"but has invalid Cell_death={row['Cell_death']!r}"
                )
            tested_counts[row["Cell_death"]] += 1

        order = plant_order[species]
        order_group = order if order in TOP_ORDERS else ORDER_OTHER_LABEL
        is_solo1 = "True" if protein in solo1_accs else "False"
        if is_solo1 == "True":
            solo1_tip_count += 1

        all_orders.add(order)
        order_group_counts[order_group] += 1

        records.append(
            {
                "tip_label": label,
                "species": species,
                "assembly": assembly,
                "protein_ID": protein,
                "plant_order": order,
                "plant_order_group": order_group,
                "sensor_clade": row["sensor_clade"],
                "is_solo1_clade": is_solo1,
                "weighted_rank": row["weighted_rank"],
                "ipTM_LCB": row["ipTM_LCB"],
                "D_APEX": row["D_APEX"],
                "screening_stage": row["screening_stage"],
                "Cell_death": row["Cell_death"],
                "WT_domesticated_autoactivity": row["WT_domesticated_autoactivity"],
                "screening_class": row["screening_class"],
            }
        )

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(columns)]
    for rec in records:
        lines.append("\t".join(rec[c] for c in columns))
    args.output_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")

    group_summary = ", ".join(
        f"{k}={order_group_counts[k]}" for k in (*TOP_ORDERS, ORDER_OTHER_LABEL)
    )
    print(
        f"Step 5.6 build_tip_annotations completed: tips={len(records)}, "
        f"unique_orders={len(all_orders)} ({group_summary}), "
        f"solo1_clade={solo1_tip_count}, "
        f"structural_ranking=Yes:{tested_counts['Yes']}/No:{tested_counts['No']}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
