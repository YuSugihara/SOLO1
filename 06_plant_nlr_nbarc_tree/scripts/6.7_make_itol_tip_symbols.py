#!/usr/bin/env python3
"""Step 6.7 script: colour the tree tips by N-terminal domain.

The architecture is read from the last underscore-delimited field of each
tip label. A tip is coloured when exactly one of the CC, TIR and RPW8 codes
is present; a CC architecture that also contains an integrated domain gets a
colour of its own. Tips whose architecture is absent or ambiguous are left
uncoloured, and the report records why.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


TIP_PATTERN = re.compile(r"(?<=[(,])([^():;,\s]+)(?=[:),;])")
ARCHITECTURE_PATTERN = re.compile(r"[A-Z]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tree", required=True, type=Path)
    parser.add_argument("--output-tip-symbols", required=True, type=Path)
    parser.add_argument("--output-report-tsv", required=True, type=Path)
    parser.add_argument("--color-rpw8", default="#DC373A")
    parser.add_argument("--color-tir", default="#84C4EF")
    parser.add_argument("--color-cc", default="#FED568")
    parser.add_argument(
        "--color-cc-with-o",
        "--color-cnlo",
        dest="color_cc_with_o",
        default="#2CA25F",
        help=(
            "Color for unambiguous CC architectures containing O "
            "(--color-cnlo is retained as a deprecated alias)"
        ),
    )
    parser.add_argument("--symbol-size", default=4, type=int)
    return parser.parse_args()


def extract_tip_ids_from_newick(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Input tree is empty: {path}")
    tips = TIP_PATTERN.findall(text)
    if not tips:
        raise ValueError(f"No tip IDs detected in tree: {path}")
    if len(tips) != len(set(tips)):
        raise ValueError("Duplicate tip IDs were detected in the input tree")
    return tips


def extract_architecture(tip_id: str) -> str | None:
    """Extract the final field under the NBARC1_architecture ID convention."""
    fields = tip_id.rsplit("_", 1)
    if len(fields) < 2:
        return None

    candidate = fields[-1].upper()

    # Simplified NLRtracker architectures contain N and/or L. Requiring one
    # prevents unannotated gene names from being mistaken for architectures.
    if not ARCHITECTURE_PATTERN.fullmatch(candidate):
        return None
    if "N" not in candidate and "L" not in candidate:
        return None
    return candidate


def classify_architecture(architecture: str | None) -> tuple[str | None, str]:
    if architecture is None:
        return None, "architecture_not_detected"

    detected = [code for code in ("C", "T", "R") if code in architecture]
    if len(detected) == 0:
        return None, "no_C_T_or_R"
    if len(detected) > 1:
        return None, f"ambiguous_{''.join(detected)}"

    class_map = {"C": "CC", "T": "TIR", "R": "RPW8"}
    return class_map[detected[0]], f"single_{detected[0]}_domain"


def write_tip_symbols(
    path: Path,
    colored_rows: list[tuple[str, str, str]],
    color_map: dict[str, str],
    symbol_size: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_SYMBOL\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tNLR_N-terminal_domain\n")
        handle.write("COLOR\t#000000\n")
        handle.write("MAXIMUM_SIZE\t12\n")
        handle.write("LEGEND_TITLE\tN-terminal domain\n")
        handle.write("LEGEND_SHAPES\t2\t2\t2\t2\n")
        handle.write(
            f"LEGEND_COLORS\t{color_map['RPW8']}\t"
            f"{color_map['TIR']}\t{color_map['CC']}\t"
            f"{color_map['CC_WITH_O']}\n"
        )
        handle.write(
            "LEGEND_LABELS\tRPW8\tTIR\tCC\tCC (architecture contains O)\n"
        )
        handle.write("DATA\n")
        for tip_id, ntd_label, color_category in colored_rows:
            symbol_label = (
                "CC (architecture contains O)"
                if color_category == "CC_WITH_O"
                else ntd_label
            )
            handle.write(
                f"{tip_id}\t2\t{symbol_size}\t"
                f"{color_map[color_category]}\t1\t1\t{symbol_label}\n"
            )


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tip_id",
        "domain_architecture",
        "ntd_label",
        "color_category",
        "colored",
        "classification_reason",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if args.symbol_size < 1:
        raise ValueError("--symbol-size must be at least 1")

    color_map = {
        "RPW8": args.color_rpw8,
        "TIR": args.color_tir,
        "CC": args.color_cc,
        "CC_WITH_O": args.color_cc_with_o,
    }
    colored_rows: list[tuple[str, str, str]] = []
    report_rows: list[dict[str, str]] = []
    class_counts: Counter[str] = Counter()

    for tip_id in extract_tip_ids_from_newick(args.input_tree):
        architecture = extract_architecture(tip_id)
        ntd_label, reason = classify_architecture(architecture)
        if ntd_label is None:
            report_rows.append(
                {
                    "tip_id": tip_id,
                    "domain_architecture": architecture or "NA",
                    "ntd_label": "NA",
                    "color_category": "NA",
                    "colored": "no",
                    "classification_reason": reason,
                }
            )
            class_counts["uncolored"] += 1
            continue

        color_category = (
            "CC_WITH_O"
            if ntd_label == "CC" and "O" in architecture
            else ntd_label
        )
        colored_rows.append((tip_id, ntd_label, color_category))
        report_rows.append(
            {
                "tip_id": tip_id,
                "domain_architecture": architecture,
                "ntd_label": ntd_label,
                "color_category": color_category,
                "colored": "yes",
                "classification_reason": reason,
            }
        )
        class_counts[ntd_label] += 1
        if color_category == "CC_WITH_O":
            class_counts["CC_with_O_color"] += 1

    write_tip_symbols(
        args.output_tip_symbols,
        colored_rows,
        color_map,
        args.symbol_size,
    )
    write_report(args.output_report_tsv, report_rows)
    print(
        "Step 6.7 iTOL annotation completed: "
        f"tips_total={len(report_rows)}, "
        f"CC={class_counts['CC']}, "
        f"TIR={class_counts['TIR']}, "
        f"RPW8={class_counts['RPW8']}, "
        f"CC_with_O_color={class_counts['CC_with_O_color']}, "
        f"uncolored={class_counts['uncolored']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
