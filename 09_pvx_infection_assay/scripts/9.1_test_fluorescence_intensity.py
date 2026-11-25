#!/usr/bin/env python3
"""Step 9.1 script: compare fluorescence intensities between the constructs.

Every construct was infiltrated into a separate site on the same leaves, so each
pairwise comparison is made within leaves with a two-sided Wilcoxon signed-rank
test on the recorded fluorescence. No transformation or normalisation is
applied. Each channel is tested on its own: the 21 pairwise comparisons of the
seven constructs are adjusted together with the Holm method, and a compact
letter display is derived from the adjusted P values at the chosen level.

The output carries the median of each construct, the median of the paired
differences and the number of leaves on which the second construct gave the
higher value, so that the direction and the size of the difference can be read
next to the P value.
"""

import argparse
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import rankdata, wilcoxon
from statsmodels.stats.multitest import multipletests

COLUMNS = ["Channel", "Receptor", "Leaf", "Mean intensity per area"]
TEST_NAME = "two-sided Wilcoxon signed-rank (paired within leaves)"
CORRECTION = "holm"


def fail(message):
    sys.exit(f"error: {message}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-tsv", required=True, type=Path,
                   help="Fluorescence in long format, one infiltration site per row.")
    p.add_argument("--output-summary-tsv", required=True, type=Path,
                   help="Where to write one row per construct and channel, with "
                        "the descriptive statistics and the letter display.")
    p.add_argument("--output-pairwise-tsv", required=True, type=Path,
                   help="Where to write one row per pairwise comparison.")
    p.add_argument("--alpha", type=float, default=0.05,
                   help="Family-wise level at which the letters are drawn "
                        "[default: %(default)s].")
    return p.parse_args()


def load(path):
    data = pd.read_csv(path, sep="\t")
    missing = [c for c in COLUMNS if c not in data.columns]
    if missing:
        fail(f"{path.name} is missing column(s): {', '.join(missing)}")
    data["Mean intensity per area"] = pd.to_numeric(
        data["Mean intensity per area"], errors="coerce")
    if data["Mean intensity per area"].isna().any():
        fail(f"{path.name} holds a value that is not a number")
    return data


def wide(data, channel):
    """One row per leaf, one column per construct, in the order of the input."""
    block = data[data["Channel"] == channel]
    order = list(block["Receptor"].drop_duplicates())
    table = block.pivot_table(index="Leaf", columns="Receptor",
                              values="Mean intensity per area", aggfunc="first")
    table = table.reindex(columns=order)
    counts = block.groupby("Receptor")["Leaf"].nunique()
    if counts.nunique() != 1:
        fail(f"{channel}: the constructs do not share the same set of leaves")
    if table.isna().any().any():
        fail(f"{channel}: a construct is missing from at least one leaf")
    leaf_order = list(block["Leaf"].drop_duplicates())
    return table.reindex(index=leaf_order)


def compare(table, channel):
    rows = []
    for first, second in combinations(table.columns, 2):
        a = table[first].to_numpy()
        b = table[second].to_numpy()
        differences = b - a
        result = wilcoxon(b, a, alternative="two-sided")
        # scipy reports min(W+, W-); report W+, the sum of the positive ranks
        nonzero = differences[differences != 0]
        ranks = rankdata(abs(nonzero))
        rows.append({
            "Channel": channel,
            "Group_1": first,
            "Group_2": second,
            "Test": TEST_NAME,
            "n_Pairs": len(differences),
            "median_1": float(pd.Series(a).median()),
            "median_2": float(pd.Series(b).median()),
            "median_of_differences": float(pd.Series(differences).median()),
            "n_Group_2_greater": int((differences > 0).sum()),
            "n_Group_2_lower": int((differences < 0).sum()),
            "n_ties": int((differences == 0).sum()),
            "W_plus": float(ranks[nonzero > 0].sum()),
            "p_value": float(result.pvalue),
        })
    return pd.DataFrame(rows)


def letters(constructs, significant, rank_by, alpha):
    """Insert and absorb: split every set that holds a significant pair."""
    sets = [set(constructs)]
    for first, second in significant:
        split = []
        for group in sets:
            if first in group and second in group:
                split.extend([group - {first}, group - {second}])
            else:
                split.append(group)
        kept = [g for g in split if g and not any(g < other for other in split)]
        sets = []
        for group in kept:
            if group not in sets:
                sets.append(group)
    sets.sort(key=lambda g: min(rank_by[c] for c in g))
    assigned = {c: "" for c in constructs}
    for index, group in enumerate(sets):
        for construct in group:
            assigned[construct] += chr(ord("a") + index)
    return {c: "".join(sorted(v)) for c, v in assigned.items()}


def summarise(table, pairwise, channel, alpha):
    medians = {c: float(table[c].median()) for c in table.columns}
    significant = [(r.Group_1, r.Group_2) for r in pairwise.itertuples()
                   if r.p_value_Holm < alpha]
    display = letters(list(table.columns), significant, medians, alpha)
    return pd.DataFrame([{
        "Channel": channel,
        "Receptor": construct,
        "n": int(table[construct].count()),
        "Mean": float(table[construct].mean()),
        "SD": float(table[construct].std(ddof=1)),
        "Median": medians[construct],
        "Minimum": float(table[construct].min()),
        "Maximum": float(table[construct].max()),
        "CLD": display[construct],
    } for construct in table.columns])


def main():
    args = parse_args()
    data = load(args.input_tsv)
    summaries, comparisons = [], []
    for channel in data["Channel"].drop_duplicates():
        table = wide(data, channel)
        pairwise = compare(table, channel)
        pairwise["p_value_Holm"] = multipletests(
            pairwise["p_value"], method=CORRECTION)[1]
        pairwise["Significant"] = pairwise["p_value_Holm"] < args.alpha
        comparisons.append(pairwise)
        summaries.append(summarise(table, pairwise, channel, args.alpha))
    summary = pd.concat(summaries, ignore_index=True)
    pairwise = pd.concat(comparisons, ignore_index=True)

    for path, frame in ((args.output_summary_tsv, summary),
                        (args.output_pairwise_tsv, pairwise)):
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False)

    print(f"Step 9.1 test_fluorescence_intensity completed: "
          f"sites={len(data)}, channels={summary['Channel'].nunique()}, "
          f"leaves={data['Leaf'].nunique()}, "
          f"comparisons={len(pairwise)}, correction={CORRECTION}, "
          f"alpha={args.alpha}, "
          f"significant={int(pairwise['Significant'].sum())}", file=sys.stderr)


if __name__ == "__main__":
    main()
