#!/usr/bin/env python3
"""Step 5.6 script: write the iTOL annotation tracks for the candidate tree.

Four tracks are written, each keyed on the tip labels of the tree:

  sensor clades   a colour strip marking the members of the sensor clades;
  screen outcome  a colour strip marking the receptors tested in the functional
                  screen;
  apex distance   a bar track of the apex distance;
  interface score a two-field bar track of the interface confidence, split at
                  0.5.

The tracks are uploaded onto the midpoint-rooted tree in iTOL.
"""

import argparse
import sys
from pathlib import Path

CLADE_COLORS = {"Pia-2/MIC1": "#2A8C82", "NRC-S": "#6B5AA6", "Pii-2": "#D69F2D"}
PHENOTYPE_COLORS = {"DV-dependent cell death": "#E4775E",
                    "Autoactive": "#83365F",
                    "No cell death": "#2C7FA3"}
APEX_COLOR = "#C1282D"
IPTM_COLORS = ("#15415C", "#B9C3CA")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tree", required=True, type=Path,
                   help="Midpoint-rooted tree whose tip labels the tracks key on.")
    p.add_argument("--table", required=True, type=Path,
                   help="Ranked and annotated table from Step 4.1.")
    p.add_argument("--output-dir", required=True, type=Path,
                   help="Directory to write the annotation files into.")
    p.add_argument("--output-prefix", default="5.6_itol",
                   help="Basename shared by the annotation files "
                        "(default: %(default)s).")
    p.add_argument("--iptm-threshold", type=float, default=0.5,
                   help="Value at which the interface-score track is split "
                        "(default: %(default)s).")
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
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        needed = ("protein_ID", "sensor_clade", "screening_class", "D_APEX", "ipTM_LCB")
        for column in needed:
            if column not in header:
                fail(f"{path} has no '{column}' column")
        index = {c: header.index(c) for c in needed}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            rows[fields[index["protein_ID"]]] = {c: fields[i] for c, i in index.items()}
    return rows


def write_file(path, header_lines, data_lines):
    path.write_text("\n".join([*header_lines, "DATA", *data_lines]) + "\n",
                    encoding="utf-8")
    return len(data_lines)


def main():
    args = parse_args()
    labels = read_tip_labels(args.tree)
    table = read_table(args.table)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Tip labels read <species>|<assembly>|<protein>; the table is keyed on the
    # protein accession alone.
    tips = {}
    for label in labels:
        parts = label.split("|")
        if len(parts) != 3:
            fail(f"tip label {label!r} is not <species>|<assembly>|<protein>")
        protein = parts[2]
        if protein not in table:
            fail(f"{protein} is a tip of the tree but is absent from {args.table}")
        tips[label] = table[protein]

    prefix = args.output_prefix
    written = {}

    clade_rows = [f"{label}\t{CLADE_COLORS[row['sensor_clade']]}\t"
                  f"{row['sensor_clade']}"
                  for label, row in tips.items() if row["sensor_clade"] in CLADE_COLORS]
    written["sensor_clades"] = write_file(
        args.output_dir / f"{prefix}_sensor_clades.txt",
        ["DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tSensor clades",
         "COLOR\t#777777", "COLOR_BRANCHES\t0", "STRIP_WIDTH\t20", "MARGIN\t2",
         "BORDER_WIDTH\t0", "SHOW_STRIP_LABELS\t0",
         "LEGEND_TITLE\tSensor clade",
         "LEGEND_SHAPES\t" + "\t".join("1" for _ in CLADE_COLORS),
         "LEGEND_COLORS\t" + "\t".join(CLADE_COLORS[k] for k in CLADE_COLORS),
         "LEGEND_LABELS\t" + "\t".join(CLADE_COLORS)],
        clade_rows)

    phenotype_rows = [f"{label}\t{PHENOTYPE_COLORS[row['screening_class']]}\t"
                      f"{row['screening_class']}"
                      for label, row in tips.items()
                      if row["screening_class"] in PHENOTYPE_COLORS]
    written["screen_outcome"] = write_file(
        args.output_dir / f"{prefix}_screen_outcome.txt",
        ["DATASET_COLORSTRIP", "SEPARATOR TAB",
         "DATASET_LABEL\tFunctional screen outcome",
         "COLOR\t#777777", "COLOR_BRANCHES\t0", "STRIP_WIDTH\t20", "MARGIN\t2",
         "BORDER_WIDTH\t0", "SHOW_STRIP_LABELS\t0", "SHOW_LABELS\t1",
         "LEGEND_TITLE\tScreen outcome",
         "LEGEND_SHAPES\t" + "\t".join("1" for _ in PHENOTYPE_COLORS),
         "LEGEND_COLORS\t" + "\t".join(PHENOTYPE_COLORS.values()),
         "LEGEND_LABELS\t" + "\t".join(PHENOTYPE_COLORS)],
        phenotype_rows)

    written["apex_distance"] = write_file(
        args.output_dir / f"{prefix}_apex_distance.txt",
        ["DATASET_SIMPLEBAR", "SEPARATOR TAB", "DATASET_LABEL\tD_APEX",
         f"COLOR\t{APEX_COLOR}", "WIDTH\t220", "MARGIN\t4", "HEIGHT_FACTOR\t0.85",
         "BAR_ZERO\t0", "SHOW_VALUE\t0", "SHOW_LABELS\t0", "BORDER_WIDTH\t0"],
        [f"{label}\t{row['D_APEX']}" for label, row in tips.items()])

    threshold = args.iptm_threshold
    iptm_rows = []
    for label, row in tips.items():
        value = float(row["ipTM_LCB"])
        above = f"{value:.10g}" if value >= threshold else "0"
        below = "0" if value >= threshold else f"{value:.10g}"
        iptm_rows.append(f"{label}\t{above}\t{below}")
    written["interface_score"] = write_file(
        args.output_dir / f"{prefix}_interface_score.txt",
        ["DATASET_MULTIBAR", "SEPARATOR TAB",
         f"DATASET_LABEL\tipTM_LCB (split at {threshold:g})",
         f"COLOR\t{IPTM_COLORS[0]}",
         "FIELD_COLORS\t" + "\t".join(IPTM_COLORS),
         f"FIELD_LABELS\tipTM >= {threshold:g}\tipTM < {threshold:g}",
         f"DATASET_SCALE\t{threshold:g}-{threshold:g}-#555555-1-1-1",
         "LEGEND_TITLE\tipTM_LCB",
         "LEGEND_SHAPES\t1\t1",
         "LEGEND_COLORS\t" + "\t".join(IPTM_COLORS),
         f"LEGEND_LABELS\t>= {threshold:g}\t< {threshold:g}",
         "WIDTH\t100", "MARGIN\t4", "HEIGHT_FACTOR\t0.85", "BORDER_WIDTH\t0",
         "SHOW_VALUE\t0", "ALIGN_FIELDS\t1"],
        iptm_rows)

    print("Step 5.6 make_itol_annotations completed: tips={}, {}".format(
        len(tips), ", ".join(f"{k}={v}" for k, v in written.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
