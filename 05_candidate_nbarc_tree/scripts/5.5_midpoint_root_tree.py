#!/usr/bin/env python3
"""Step 5.5 script: place the root of the tree at its midpoint.

Midpoint-root the maximum-likelihood tree with Biopython and write it with
branch lengths to five decimals.
"""

import argparse
import sys
from pathlib import Path

from Bio import Phylo

BRANCH_LENGTH_FORMAT = "%.5f"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-tree", required=True, type=Path,
                   help="Newick tree to root.")
    p.add_argument("--output-tree", required=True, type=Path,
                   help="Midpoint-rooted Newick tree.")
    return p.parse_args()


def main():
    args = parse_args()
    tree = Phylo.read(str(args.input_tree), "newick")
    n_tips = len(tree.get_terminals())
    if n_tips < 3:
        sys.exit(f"error: {args.input_tree} has {n_tips} tips; rooting needs at least 3")
    tree.root_at_midpoint()
    if len(tree.get_terminals()) != n_tips:
        sys.exit("error: rooting changed the number of tips")
    args.output_tree.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tree.open("w", encoding="utf-8") as handle:
        Phylo.write(tree, handle, "newick",
                    format_branch_length=BRANCH_LENGTH_FORMAT)
    print(f"Step 5.5 midpoint_root_tree completed: tips={n_tips}", file=sys.stderr)


if __name__ == "__main__":
    main()
