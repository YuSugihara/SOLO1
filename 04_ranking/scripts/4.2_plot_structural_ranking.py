#!/usr/bin/env python3
"""Step 4.2 script: draw the structural ranking panel.

Entries are placed along the x axis in the order of the weighted ranking, with
the rank scale along the top. The apex distance and the interface confidence
are drawn as bar tracks, and one lane per sensor clade marks where the members
of each clade fall. Representative ranks are marked with dashed lines.
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BAR_COLORS = {"D_APEX": "#C45B47", "ipTM_LCB": "#15415C"}
CLADE_LANES = [("Pia-2/MIC1", "Pia-2/MIC1", "#2A8C82"),
               ("NRC-S", "NRC-S", "#6B5AA6"),
               ("Pii-2", "Pii-2", "#D69F2D")]
GUIDE_COLOR = "#374151"
LANE_BASELINE = "#D2D6DC"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True, type=Path,
                   help="Ranked and annotated table from Step 4.1.")
    p.add_argument("--output-prefix", required=True, type=Path,
                   help="Path without extension; PDF, PNG and SVG are written.")
    p.add_argument("--reference-id", default="NbNRC2a",
                   help="Entry marked as the benchmark (default: %(default)s).")
    p.add_argument("--highlight-ranks", nargs="+", type=int, default=[1, 18, 128, 297],
                   help="Ranks marked with vertical dashed guide lines (default: %(default)s).")
    p.add_argument("--dpi", type=int, default=600,
                   help="Resolution of the raster output (default: %(default)s).")
    return p.parse_args()


def fail(message):
    sys.exit(f"error: {message}")


def configure_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.default": "regular",
        "axes.linewidth": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def draw_bar_track(ax, ranks, values, color, label, x_max, y_scale=1.08, top_axis=False):
    ax.bar(ranks, values, width=1.0, color=color, edgecolor="none", zorder=2)
    ax.set_xlim(0.5, x_max + 0.5)
    ax.set_ylim(0, float(values.max()) * y_scale)
    ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=8, fontsize=8.5)
    ax.set_yticks([])
    for side in ("right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    if top_axis:
        ax.spines["top"].set_visible(True)
        ax.spines["top"].set_linewidth(0.75)
        ax.spines["top"].set_color("#1A1A1A")
        ax.xaxis.set_ticks_position("top")
    else:
        ax.spines["top"].set_visible(False)
        ax.tick_params(axis="x", which="both", top=False, bottom=False,
                       labeltop=False, labelbottom=False, length=0)


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
    ax.set_ylim(len(CLADE_LANES) - 0.45, -0.82)
    ax.set_yticks([])
    ax.tick_params(axis="x", which="both", top=False, bottom=False,
                   labeltop=False, labelbottom=False, length=0)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.text(x_max + 0.5 + x_max * 0.012, -0.76, "Sensor clades",
            fontsize=8, va="center", ha="left", fontweight="bold", color="#1A1A1A")


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

    highlight_ranks = sorted(set([r for r in args.highlight_ranks if 1 <= r <= x_max]))

    fig = plt.figure(figsize=(8.6, 1.85), facecolor="white")
    grid = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[0.60, 0.54, 1.25],
                            hspace=0.07, left=0.085, right=0.845,
                            top=0.72, bottom=0.05)
    ax_apex = fig.add_subplot(grid[0])
    ax_iptm = fig.add_subplot(grid[1], sharex=ax_apex)
    ax_clade = fig.add_subplot(grid[2], sharex=ax_apex)

    ranks = data["weighted_rank"].to_numpy()
    draw_bar_track(ax_apex, ranks, data["D_APEX"], BAR_COLORS["D_APEX"],
                   r"$\mathrm{D_{APEX}}$ (Å)", x_max, y_scale=1.12, top_axis=True)
    draw_bar_track(ax_iptm, ranks, data["ipTM_LCB"], BAR_COLORS["ipTM_LCB"],
                   "ipTM", x_max, y_scale=1.06, top_axis=False)
    draw_clade_lanes(ax_clade, data, x_max)

    # Mark the representative ranks across all three tracks.
    for r in highlight_ranks:
        for ax in (ax_apex, ax_iptm, ax_clade):
            ax.axvline(r, color=GUIDE_COLOR, lw=0.75, ls=(0, (2.6, 2.0)),
                       alpha=0.80, zorder=3, clip_on=False)

    # Top rank axis ticks and label.
    base_ticks = [t for t in (50, 100, 150, 200, 250) if t <= x_max]
    ticks = sorted(set([1, reference_rank, x_max] + highlight_ranks + base_ticks))
    ax_apex.set_xticks(ticks)
    ax_apex.set_xticklabels([str(t) for t in ticks], fontsize=7.8)
    ax_apex.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False,
                        length=3.2, width=0.75, pad=2.2, colors="#4B5563")
    for lbl, tval in zip(ax_apex.get_xticklabels(), ticks):
        if tval in highlight_ranks:
            lbl.set_fontweight("bold")
            lbl.set_color("#111827")
            lbl.set_fontsize(8.0)

    ax_apex.xaxis.set_label_position("top")
    ax_apex.set_xlabel("Composite structural rank", fontsize=9, labelpad=7, color="#1A1A1A")

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
