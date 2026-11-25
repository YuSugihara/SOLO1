#!/usr/bin/env python3
"""Step 4.3 script: plot the predicted aligned error of an AlphaFold 3 multimer prediction.

The matrix covers every token of the prediction, ligands included, so the
confidence in the placement of the membrane proxy is visible alongside the
protein chains.

Borders are drawn between protein chains and once where the protein chains
end, never between individual ligand copies: a prediction carrying tens of
ligand copies would otherwise be covered by border lines. Each protein chain
carries its number above the matrix and the ligand block carries a word label,
so the panel reads without residue numbers or axis titles; both can be switched
on.

Protein chains are told apart from ligand copies by their tokens: a protein
chain spans many residue identifiers, whereas each ligand copy is a single
residue whose atoms are tokenised individually.
"""

import argparse
import json
import sys
from itertools import groupby
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

PAE_CMAP = LinearSegmentedColormap.from_list(
    "af3_reference_green",
    ["#087f3a", "#2fa24f", "#74c476", "#c7e9c0", "#f7fcf5"],
    N=256,
)
FONT_STACK = ["Helvetica", "Arial", "DejaVu Sans"]
INK = "#222222"


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "full_data_json", type=Path,
        help="AlphaFold 3 full_data JSON holding 'pae', 'token_chain_ids' and 'token_res_ids'.",
    )
    parser.add_argument("--output-prefix", default=None,
        help="Output prefix. Defaults to '<full_data_json stem>_pae'.")
    parser.add_argument("--vmax", type=float, default=31.0,
        help="Upper bound of the colour scale, in angstroms (default %(default)s).")
    parser.add_argument("--ligand-label", default="oleic acids",
        help="Word label written above the ligand block (default '%(default)s').")
    parser.add_argument("--ticks", action="store_true",
        help="Write token numbers on both axes. Omitted by default.")
    parser.add_argument("--axis-titles", action="store_true",
        help="Write the axis titles. Omitted by default.")
    parser.add_argument("--colorbar", action="store_true",
        help="Draw the colour bar beneath the matrix. Omitted by default.")
    parser.add_argument("--no-chain-borders", action="store_true",
        help="Leave out the borders between protein chains.")
    parser.add_argument("--no-chain-labels", action="store_true",
        help="Leave out the chain numbers and the ligand label above the matrix.")
    parser.add_argument("--dpi", type=int, default=200,
        help="Resolution of raster outputs (default %(default)s).")
    parser.add_argument("--formats", nargs="+", default=["png"],
        help="Output formats to write (default %(default)s).")
    return parser.parse_args()


def read_pae(path):
    with path.open() as handle:
        data = json.load(handle)
    for key in ("pae", "token_chain_ids", "token_res_ids"):
        if key not in data:
            raise SystemExit(f"{path}: '{key}' is missing from the JSON")
    pae = np.asarray(data["pae"], dtype=np.float32)
    chain_ids = list(data["token_chain_ids"])
    res_ids = list(data["token_res_ids"])
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise SystemExit(f"{path}: the PAE matrix is not square; got {pae.shape}")
    if len(chain_ids) != pae.shape[0] or len(res_ids) != pae.shape[0]:
        raise SystemExit(
            f"{path}: token count does not match the matrix: {len(chain_ids)} chain ids, "
            f"{len(res_ids)} residue ids, matrix {pae.shape[0]}"
        )
    return pae, chain_ids, res_ids


def split_chains(chain_ids, res_ids):
    """Return (chain_id, start, end, is_protein) per chain, in token order."""
    chains, start = [], 0
    for chain_id, group in groupby(chain_ids):
        end = start + len(list(group))
        chains.append((chain_id, start, end, len(set(res_ids[start:end])) > 1))
        start = end
    return chains


def annotate_blocks(ax, chains, n_tokens, ligand_label):
    """Number each protein chain above the matrix and label the ligand block."""
    baseline = -0.5 - n_tokens * 0.012
    protein = [c for c in chains if c[3]]
    for number, (_, start, end, _) in enumerate(protein, start=1):
        ax.text((start + end - 1) / 2, baseline, str(number), ha="center", va="bottom",
                fontsize=8, color=INK, clip_on=False)
    ligand = [c for c in chains if not c[3]]
    if ligand:
        start, end = ligand[0][1], ligand[-1][2]
        ax.text((start + end - 1) / 2, baseline, ligand_label, ha="center", va="bottom",
                fontsize=8, color=INK, clip_on=False)


def plot(pae, chains, prefix, args):
    n_tokens = pae.shape[0]
    protein = [c for c in chains if c[3]]
    borders = ([end for _, _, end, is_protein in chains
                if is_protein and end < n_tokens] if not args.no_chain_borders else [])

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
        "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })

    if args.colorbar:
        fig = plt.figure(figsize=(4.9, 5.6), facecolor="white")
        grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 0.055], hspace=0.30)
        ax = fig.add_subplot(grid[0])
    else:
        fig = plt.figure(figsize=(4.9, 5.0), facecolor="white")
        ax = fig.add_subplot(1, 1, 1)

    image = ax.imshow(pae, origin="upper", cmap=PAE_CMAP, vmin=0, vmax=args.vmax,
                      interpolation="nearest", aspect="equal")
    for position in borders:
        ax.axhline(position - 0.5, color=INK, linewidth=0.65)
        ax.axvline(position - 0.5, color=INK, linewidth=0.65)

    ax.set_xlim(-0.5, n_tokens - 0.5)
    ax.set_ylim(n_tokens - 0.5, -0.5)

    if args.ticks:
        positions = [0] + [end - 1 for _, _, end, is_protein in chains if is_protein]
        if positions[-1] != n_tokens - 1:
            positions.append(n_tokens - 1)
        ax.set_xticks(positions); ax.set_yticks(positions)
        ax.set_xticklabels([str(p + 1) for p in positions])
        ax.set_yticklabels([str(p + 1) for p in positions])
        ax.tick_params(length=0, pad=6, colors=INK)
    else:
        ax.set_xticks([]); ax.set_yticks([])

    if args.axis_titles:
        ax.set_xlabel("Scored token", labelpad=12)
        ax.set_ylabel("Aligned token", labelpad=10)

    for spine in ax.spines.values():
        spine.set_visible(True); spine.set_color(INK); spine.set_linewidth(0.65)

    if not args.no_chain_labels:
        annotate_blocks(ax, chains, n_tokens, args.ligand_label)

    if args.colorbar:
        cax = fig.add_subplot(grid[1])
        bar = fig.colorbar(image, cax=cax, orientation="horizontal")
        bar.set_ticks(np.arange(0, int(args.vmax) + 1, 5))
        bar.set_label("Expected Position Error (Å)", labelpad=4)
        bar.outline.set_linewidth(0.6)
        bar.ax.tick_params(axis="x", direction="in", length=4, width=0.7, pad=6, colors=INK)

    written = []
    for file_format in args.formats:
        path = Path(f"{prefix}.{file_format}")
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight", facecolor="white")
        written.append(path)
    plt.close(fig)
    return written, protein


def main():
    args = parse_args()
    prefix = args.output_prefix or str(
        args.full_data_json.with_name(f"{args.full_data_json.stem}_pae"))
    pae, chain_ids, res_ids = read_pae(args.full_data_json)
    chains = split_chains(chain_ids, res_ids)
    written, protein = plot(pae, chains, prefix, args)
    for path in written:
        print(f"Saved: {path}")
    print(
        f"Step 4.3 plot_af3_pae_multimer completed: tokens={pae.shape[0]}, "
        f"protein_chains={len(protein)}, ligand_copies={len(chains) - len(protein)}, "
        f"pae_max={float(np.nanmax(pae)):.2f}, "
        f"borders={'no' if args.no_chain_borders else 'yes'}, "
        f"labels={'no' if args.no_chain_labels else 'yes'}, "
        f"ticks={'yes' if args.ticks else 'no'}, colorbar={'yes' if args.colorbar else 'no'}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
