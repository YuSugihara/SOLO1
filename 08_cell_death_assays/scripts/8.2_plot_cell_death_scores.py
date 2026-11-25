#!/usr/bin/env python3
"""Step 8.2 script: draw the cell death score panel of an assay.

Each receptor occupies one block on a shared y axis of scores. Within a block
the two conditions are drawn side by side, and within a condition each
experiment date is drawn as its own column. Observations that share a score
within a column are collapsed into one point whose area is proportional to the
number of observations.

The P value of each receptor is read from the table written by Step 8.1 and is
drawn above its block; the script performs no test of its own.
"""

import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SCORE_COLUMNS = ["NLR", "Experimental_Condition", "Experiment_Date", "Cell death score"]
SCORE_MIN, SCORE_MAX = 0, 7
POINT_AREA = 30.0
REPLICATE_SPACING, CONDITION_GAP, RECEPTOR_GAP = 0.52, 0.38, 0.82


def fail(message):
    sys.exit(f"error: {message}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-tsv", required=True, type=Path,
                   help="Scores in long format, as read by Step 8.1.")
    p.add_argument("--statistics-tsv", required=True, type=Path,
                   help="Test results written by Step 8.1.")
    p.add_argument("--output-prefix", required=True, type=Path,
                   help="Path without extension; PDF, SVG and PNG are written.")
    p.add_argument("--colours", nargs=2, metavar=("FIRST", "SECOND"),
                   default=("#3070AD", "#C66526"),
                   help="One colour per condition (default: %(default)s).")
    p.add_argument("--width-scale", type=float, default=1.0,
                   help="Multiplier on the calculated figure width (default: %(default)s).")
    p.add_argument("--dpi", type=int, default=600,
                   help="Resolution of the raster output (default: %(default)s).")
    return p.parse_args()


def configure_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def p_value_label(p_value):
    for threshold in (0.001, 0.01, 0.05):
        if p_value < threshold:
            return f"$p < {threshold:g}$"
    return "ns"


def load(scores_path, statistics_path):
    if not scores_path.is_file():
        fail(f"input not found: {scores_path}")
    if not statistics_path.is_file():
        fail(f"statistics not found: {statistics_path}")
    data = pd.read_csv(scores_path, sep="\t")
    missing = [c for c in SCORE_COLUMNS if c not in data.columns]
    if missing:
        fail(f"{scores_path.name} is missing column(s): {', '.join(missing)}")
    data = data.loc[:, SCORE_COLUMNS].dropna(how="all")
    data["Experiment_Date"] = data["Experiment_Date"].astype("string").str.strip()
    data["Cell death score"] = pd.to_numeric(data["Cell death score"], errors="coerce")
    if data.isna().any(axis=None):
        fail(f"{scores_path.name} has missing or unreadable values")
    scores = data["Cell death score"]
    if scores.lt(SCORE_MIN).any() or scores.gt(SCORE_MAX).any():
        fail(f"scores must lie on the {SCORE_MIN} to {SCORE_MAX} scale")

    statistics = pd.read_csv(statistics_path, sep="\t")
    for column in ("NLR", "Experimental_Condition_1", "Experimental_Condition_2", "p_value"):
        if column not in statistics.columns:
            fail(f"{statistics_path.name} is missing column {column}")
    conditions = (str(statistics["Experimental_Condition_1"].iloc[0]),
                  str(statistics["Experimental_Condition_2"].iloc[0]))
    if statistics["Experimental_Condition_1"].nunique() != 1 or \
       statistics["Experimental_Condition_2"].nunique() != 1:
        fail(f"{statistics_path.name} mixes more than one pair of conditions")
    present = set(data["Experimental_Condition"])
    if present != set(conditions):
        fail(f"conditions differ between the two files: {sorted(present)} and {list(conditions)}")
    receptors = list(statistics["NLR"])
    if sorted(receptors) != sorted(data["NLR"].unique()):
        fail("the two files cover different receptors")
    data["NLR"] = pd.Categorical(data["NLR"], categories=receptors, ordered=True)
    return data, statistics.set_index("NLR"), conditions


def lay_out(data, receptors, dates_of):
    layout, cursor = {}, 0.0
    for receptor in receptors:
        positions, start = {}, cursor
        for index, condition in enumerate(("first", "second")):
            columns = [cursor + i * REPLICATE_SPACING for i in range(len(dates_of[receptor]))]
            positions[condition] = columns
            cursor = columns[-1] + REPLICATE_SPACING
            if index == 0:
                cursor += CONDITION_GAP
        end = positions["second"][-1]
        layout[receptor] = {"columns": positions,
                            "left": start - REPLICATE_SPACING * 0.38,
                            "right": end + REPLICATE_SPACING * 0.38,
                            "centre": (start + end) / 2}
        cursor = end + RECEPTOR_GAP
    return layout


def draw(data, statistics, conditions, colours, width_scale, dpi, output_prefix):
    receptors = [str(v) for v in data["NLR"].cat.categories]
    dates_of = {r: sorted(data.loc[data["NLR"].astype("string") == r, "Experiment_Date"].unique())
                for r in receptors}
    layout = lay_out(data, receptors, dates_of)
    colour_of = dict(zip(("first", "second"), colours))
    condition_of = dict(zip(("first", "second"), conditions))

    x_min = float(layout[receptors[0]]["left"])
    x_max = float(layout[receptors[-1]]["right"])
    width = max(7.5, 1.4 + (x_max - x_min) * 0.38) * width_scale
    fig, ax = plt.subplots(figsize=(width, 5.4))

    bracket_y, bracket_height = SCORE_MAX + 0.52, 0.20
    y_lower, y_upper = SCORE_MIN - 0.35, SCORE_MAX + 1.65
    wrap = max(12, round(28 * min(1.0, max(len(d) for d in dates_of.values()) / 4)
                         * min(width_scale, 1.0) ** 3))
    labels = []

    for index, receptor in enumerate(receptors):
        subset = data.loc[data["NLR"].astype("string") == receptor]
        block = layout[receptor]
        for slot in ("first", "second"):
            for column, date in zip(block["columns"][slot], dates_of[receptor]):
                values = subset.loc[(subset["Experiment_Date"] == date) &
                                    (subset["Experimental_Condition"] == condition_of[slot]),
                                    "Cell death score"].to_numpy()
                if len(values) == 0:
                    continue
                unique, counts = np.unique(values, return_counts=True)
                ax.scatter(np.full(len(unique), column), unique, s=POINT_AREA * counts,
                           facecolor=colour_of[slot], edgecolor="none", zorder=3)

        left = float(np.mean(block["columns"]["first"]))
        right = float(np.mean(block["columns"]["second"]))
        ax.plot([left, left, right, right],
                [bracket_y, bracket_y + bracket_height, bracket_y + bracket_height, bracket_y],
                color="black", linewidth=0.8, clip_on=False)
        ax.text((left + right) / 2, bracket_y + bracket_height + 0.08,
                p_value_label(float(statistics.loc[receptor, "p_value"])),
                ha="center", va="bottom", fontsize=9)

        labels.append(ax.text(block["centre"], -0.065, "\n".join(textwrap.wrap(receptor, wrap)),
                              transform=ax.get_xaxis_transform(), ha="center", va="top",
                              fontsize=8.5, fontweight="bold", clip_on=False))
        if index < len(receptors) - 1:
            boundary = (float(block["right"]) + float(layout[receptors[index + 1]]["left"])) / 2
            ax.axvline(boundary, color="#D0D0D0", linewidth=0.7, linestyle=(0, (2, 3)), zorder=0)

    ax.set_xlim(x_min - 0.15, x_max + 0.15)
    ax.set_ylim(y_lower, y_upper)
    ax.set_yticks(range(SCORE_MIN, SCORE_MAX + 1))
    ax.set_xticks([])
    ax.set_ylabel("Cell death score")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_bounds(y_lower, SCORE_MAX)

    condition_legend = ax.legend(
        handles=[Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=colour_of[s],
                        markeredgecolor="none", markersize=7, label=condition_of[s])
                 for s in ("first", "second")],
        loc="upper left", bbox_to_anchor=(1.002, 0.72), frameon=False, borderaxespad=0)
    ax.add_artist(condition_legend)

    largest = int(data.groupby(["NLR", "Experimental_Condition", "Experiment_Date",
                               "Cell death score"], observed=True).size().max())
    steps = sorted({1, max(1, -(-largest // 2)), largest}, reverse=True)
    count_legend = ax.legend(
        handles=[ax.scatter([], [], s=POINT_AREA * c, facecolor="black",
                            edgecolor="none", label=str(c)) for c in steps],
        title="Count", loc="upper left", bbox_to_anchor=(1.002, 0.47), frameon=False,
        borderaxespad=0, labelspacing=0.9, handletextpad=0.8)

    fig.subplots_adjust(left=0.07, right=0.93, top=0.88, bottom=0.19)
    extra = (condition_legend, count_legend, *labels)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in ((".pdf", {}), (".svg", {}), (".png", {"dpi": dpi})):
        fig.savefig(Path(str(output_prefix) + suffix), bbox_inches="tight",
                    bbox_extra_artists=extra, pad_inches=0.18, facecolor="white", **kwargs)
    plt.close(fig)
    return len(receptors)


def main():
    args = parse_args()
    configure_style()
    data, statistics, conditions = load(args.input_tsv, args.statistics_tsv)
    receptors = draw(data, statistics, conditions, tuple(args.colours),
                     args.width_scale, args.dpi, args.output_prefix)
    print(f"Step 8.2 plot_cell_death_scores completed: observations={len(data)}, "
          f"receptors={receptors}, conditions={conditions[0]} vs {conditions[1]}, "
          f"dates={data['Experiment_Date'].nunique()}", file=sys.stderr)


if __name__ == "__main__":
    main()
