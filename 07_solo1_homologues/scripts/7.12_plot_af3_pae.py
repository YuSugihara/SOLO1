#!/usr/bin/env python3
"""Step 7.12 script: plot the predicted aligned error of an AlphaFold 3 model.

Reads the predicted aligned error (PAE) matrix and the token chain assignments
from an AlphaFold 3 `full_data` JSON and draws the matrix as a heat map, with a
line at each chain boundary and a horizontal colour bar. The colour scale runs
from 0 to the value of --vmax in angstroms.

The `full_data` JSON files are part of the AlphaFold 3 output and are too large
for this repository; they are deposited separately. The job requests that
produced them are in `data/inputs/7.0_af3_job_requests`.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


PAE_CMAP = LinearSegmentedColormap.from_list(
    "af3_reference_green",
    ["#087f3a", "#2fa24f", "#74c476", "#c7e9c0", "#f7fcf5"],
    N=256,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot an AlphaFold 3 PAE matrix with chain-border lines."
    )
    parser.add_argument(
        "full_data_json",
        help="AlphaFold 3 full_data JSON containing 'pae' and 'token_chain_ids'.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix. Defaults to '<full_data_json stem>_pae_plot'.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=31.0,
        help="Maximum PAE value for the color scale.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Resolution for raster outputs.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Output formats to write.",
    )
    return parser.parse_args()


def read_pae_data(full_data_json):
    with full_data_json.open() as handle:
        data = json.load(handle)

    pae = np.asarray(data["pae"], dtype=np.float32)
    token_chain_ids = np.asarray(data["token_chain_ids"])

    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE matrix must be square; got {pae.shape}")
    if len(token_chain_ids) != pae.shape[0]:
        raise ValueError(
            "token_chain_ids length does not match PAE size: "
            f"{len(token_chain_ids)} vs {pae.shape[0]}"
        )

    return pae, token_chain_ids


def get_chain_boundaries(token_chain_ids):
    boundaries = []
    start = 0
    for index in range(1, len(token_chain_ids) + 1):
        if index == len(token_chain_ids) or token_chain_ids[index] != token_chain_ids[start]:
            boundaries.append((start, index, str(token_chain_ids[start])))
            start = index
    return boundaries


def plot_pae(pae, chain_boundaries, output_prefix, vmax=31.0, dpi=600, formats=None):
    if formats is None:
        formats = ["png", "pdf"]

    n_residues = pae.shape[0]
    tick_positions = [0] + [end - 1 for _, end, _ in chain_boundaries]
    tick_labels = [str(position + 1) for position in tick_positions]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    fig = plt.figure(figsize=(4.9, 5.85), facecolor="white")
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 0.055], hspace=0.36)
    ax = fig.add_subplot(gs[0])

    image = ax.imshow(
        pae,
        origin="upper",
        cmap=PAE_CMAP,
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal",
    )

    for _, chain_end, _ in chain_boundaries[:-1]:
        border_position = chain_end - 0.5
        ax.axhline(border_position, color="#222222", linewidth=0.65)
        ax.axvline(border_position, color="#222222", linewidth=0.65)

    ax.set_xlim(-0.5, n_residues - 0.5)
    ax.set_ylim(n_residues - 0.5, -0.5)
    ax.set_xlabel("Scored Residue", labelpad=18)
    ax.set_ylabel("Aligned Residue", labelpad=12)
    ax.set_xticks(tick_positions)
    ax.set_yticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_yticklabels(tick_labels)
    ax.tick_params(length=0, pad=7, colors="#222222")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#222222")
        spine.set_linewidth(0.65)

    cax = fig.add_subplot(gs[1])
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    cbar.set_ticks(np.arange(0, 31, 5))
    cbar.set_label("Expected Position Error (Å)", labelpad=4)
    cbar.outline.set_linewidth(0.6)
    cbar.ax.tick_params(
        axis="x",
        which="major",
        direction="in",
        length=4,
        width=0.7,
        pad=8,
        colors="#222222",
    )

    fig.subplots_adjust(left=0.16, right=0.97, top=0.98, bottom=0.10)

    output_paths = []
    for file_format in formats:
        output_path = Path(f"{output_prefix}.{file_format}")
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        output_paths.append(output_path)

    plt.close(fig)
    return output_paths


def main():
    args = parse_args()
    full_data_json = Path(args.full_data_json)
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = str(full_data_json.with_name(f"{full_data_json.stem}_pae_plot"))

    pae, token_chain_ids = read_pae_data(full_data_json)
    chain_boundaries = get_chain_boundaries(token_chain_ids)
    output_paths = plot_pae(
        pae,
        chain_boundaries,
        output_prefix,
        vmax=args.vmax,
        dpi=args.dpi,
        formats=args.formats,
    )

    print(f"PAE shape: {pae.shape}")
    print(
        "Chain boundaries: "
        + ", ".join(
            f"{chain_id}:{start + 1}-{end}"
            for start, end, chain_id in chain_boundaries
        )
    )
    print(f"PAE min/max: {float(np.nanmin(pae)):.3f}/{float(np.nanmax(pae)):.3f}")
    for output_path in output_paths:
        print(f"Saved: {output_path}")

    print(f"Step 7.12 plot_af3_pae completed: "
          f"residues={pae.shape[0]}, chains={len(chain_boundaries)}, "
          f"pae_max={float(np.nanmax(pae)):.2f}, vmax={args.vmax}, "
          f"outputs={len(output_paths)}", file=sys.stderr)


if __name__ == "__main__":
    main()
