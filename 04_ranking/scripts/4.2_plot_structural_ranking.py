#!/usr/bin/env python3
"""Step 4.2 script: draw the structural ranking panel.

Entries are placed along the x axis in the order of the weighted ranking. The
apex distance and the interface confidence are drawn as bar tracks above and
below, and one lane per sensor clade marks where the members of each clade
fall. The benchmark is marked with a dashed line and the region that ranks
above it is shaded.
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

BAR_COLORS = {"D_APEX": "#C1282D", "ipTM_LCB": "#15415C"}
CLADE_LANES = [("Pia-2/MIC1", "Pia-2/MIC1", "#2A8C82"),
               ("NRC-S", "NRC-S", "#6B5AA6"),
               ("Pii-2", "Pii-2", "#D69F2D")]
REFERENCE_COLOR = "#1A1A1A"
LANE_BASELINE = "#C8CCD0"
ABOVE_REFERENCE_FILL = "#E4EDF4"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True, type=Path,
                   help="Ranked and annotated table from Step 4.1.")
    p.add_argument("--output-prefix", required=True, type=Path,
                   help="Path without extension; PDF, PNG and SVG are written.")
    p.add_argument("--reference-id", default="NbNRC2a",
                   help="Entry marked as the benchmark (default: %(default)s).")
    p.add_argument("--dpi", type=int, default=600,
                   help="Resolution of the raster output (default: %(default)s).")
    return p.parse_args()


def fail(message):
    sys.exit(f"error: {message}")


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


def draw_bar_track(ax, ranks, values, color, label, x_max):
    ax.bar(ranks, values, width=1.0, color=color, edgecolor="none")
    ax.set_xlim(0.5, x_max + 0.5)
    ax.set_ylim(0, float(values.max()) * 1.08)
    ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=8, fontsize=8.5)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0, labelbottom=False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)


def draw_clade_lanes(ax, data, x_max):
    for y, (key, display, color) in enumerate(CLADE_LANES):
        ax.hlines(y, 0.5, x_max + 0.5, color=LANE_BASELINE, lw=0.5, zorder=1)
        selected = data.loc[data["sensor_clade"].eq(key), "weighted_rank"]
        if selected.empty:
            fail(f"no entries are annotated as the {key} clade")
        # Drawn as bars of width 1.0, matching the bars in the tracks above.
        ax.bar(selected, height=0.68, bottom=y - 0.34, width=1.0,
               color=color, edgecolor="none", zorder=2)
        ax.text(x_max + 0.5 + x_max * 0.012, y, display, color=color,
                fontsize=7.6, va="center", ha="left")
    ax.set_ylim(len(CLADE_LANES) - 0.5, -0.8)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.text(x_max + 0.5 + x_max * 0.012, -0.72, "Sensor clades",
            fontsize=8, va="center", ha="left", fontweight="bold")


def main():
    args = parse_args()
    configure_style()

    data = pd.read_csv(args.table, sep="\t", keep_default_na=False)
    for column in ("weighted_rank", "D_APEX", "ipTM_LCB", "sensor_clade", "protein_ID"):
        if column not in data.columns:
            fail(f"{args.table} has no '{column}' column")
    data = data.sort_values("weighted_rank").reset_index(drop=True)
    x_max = int(data["weighted_rank"].max())

    reference = data.loc[data["protein_ID"].eq(args.reference_id)]
    if len(reference) != 1:
        fail(f"expected exactly one entry for {args.reference_id}")
    reference_rank = int(reference["weighted_rank"].iloc[0])
    n_above = reference_rank - 1
    if n_above < 1:
        fail(f"{args.reference_id} ranks first, so nothing lies above it")

    fig = plt.figure(figsize=(8.6, 1.78), facecolor="white")
    grid = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[0.55, 0.55, 1.25],
                            hspace=0.06, left=0.085, right=0.845,
                            top=0.86, bottom=0.20)
    ax_apex = fig.add_subplot(grid[0])
    ax_iptm = fig.add_subplot(grid[1], sharex=ax_apex)
    ax_clade = fig.add_subplot(grid[2], sharex=ax_apex)

    ranks = data["weighted_rank"].to_numpy()
    draw_bar_track(ax_apex, ranks, data["D_APEX"], BAR_COLORS["D_APEX"],
                   "$D_{APEX}$ (Å)", x_max)
    draw_bar_track(ax_iptm, ranks, data["ipTM_LCB"], BAR_COLORS["ipTM_LCB"],
                   "ipTM", x_max)
    draw_clade_lanes(ax_clade, data, x_max)

    # Shade the region that ranks above the benchmark.
    for ax in (ax_apex, ax_iptm, ax_clade):
        ax.add_patch(Rectangle((0.5, 0), n_above, 1, transform=ax.get_xaxis_transform(),
                               facecolor=ABOVE_REFERENCE_FILL, edgecolor="none", zorder=0))
        ax.axvline(reference_rank, color=REFERENCE_COLOR, lw=1.0,
                   ls=(0, (3.0, 2.2)), zorder=3)

    ax_apex.plot([reference_rank], [1.0], marker="v", markersize=4.5,
                 color="black", transform=ax_apex.get_xaxis_transform(),
                 clip_on=False, zorder=4)
    ax_apex.text(reference_rank, 1.10, args.reference_id, fontsize=8,
                 ha="center", va="bottom", transform=ax_apex.get_xaxis_transform())

    # The benchmark carries a tick of its own.
    ticks = [1, reference_rank] + [t for t in (50, 100, 150, 200, 250)
                                   if t <= x_max] + [x_max]
    ax_clade.set_xticks(sorted(set(ticks)))
    ax_clade.set_xticklabels([str(t) for t in sorted(set(ticks))], fontsize=8)
    ax_clade.tick_params(axis="x", length=3, pad=2)
    ax_clade.set_xlabel("Composite structural rank", fontsize=9, labelpad=5)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    # Build the output names by appending, not with Path.with_suffix: the file
    # names carry a step number such as 4.2, which with_suffix would strip.
    for suffix, kwargs in ((".pdf", {}), (".svg", {}), (".png", {"dpi": args.dpi})):
        fig.savefig(Path(str(args.output_prefix) + suffix), facecolor="white", **kwargs)
    plt.close(fig)

    print(f"Step 4.2 plot_structural_ranking completed: entries={x_max}, "
          f"{args.reference_id}_rank={reference_rank}, above_reference={n_above}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
