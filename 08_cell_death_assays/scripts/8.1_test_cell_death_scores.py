#!/usr/bin/env python3
"""Step 8.1 script: compare cell death scores between two conditions.

Each receptor is tested on its own. The two conditions are compared with a
two-sided Mann-Whitney U test on the pooled observations, using the asymptotic
approximation with the tie correction and the continuity correction, as the
scores are an ordinal 0 to 7 scale with many ties. The experiment date is read
and reported, but it takes no part in the test: the observations of all dates
are pooled.

The output carries the sample size and the median of each condition, the U
statistic of the first condition, and the probability of superiority, so that
the size of the difference can be read next to the P value.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

COLUMNS = ["NLR", "Experimental_Condition", "Experiment_Date", "Cell death score"]
TEST_NAME = "two-sided Mann-Whitney U (asymptotic, tie-corrected)"
SCORE_MIN, SCORE_MAX = 0, 7


def fail(message):
    sys.exit(f"error: {message}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-tsv", required=True, type=Path,
                   help="Scores in long format, one observation per row.")
    p.add_argument("--output-tsv", required=True, type=Path,
                   help="Where to write one row of test results per receptor.")
    p.add_argument("--conditions", nargs=2, metavar=("FIRST", "SECOND"), default=None,
                   help="The two conditions to compare, in the order they are "
                        "reported. Defaults to their order of first appearance.")
    return p.parse_args()


def significance(p_value):
    for threshold, label in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")):
        if p_value < threshold:
            return label
    return "ns"


def load(path, requested):
    if not path.is_file():
        fail(f"input not found: {path}")
    data = pd.read_csv(path, sep="\t")
    missing = [c for c in COLUMNS if c not in data.columns]
    if missing:
        fail(f"{path.name} is missing column(s): {', '.join(missing)}")
    data = data.loc[:, COLUMNS].dropna(how="all")
    data["NLR"] = data["NLR"].astype("string").str.strip()
    data["Experimental_Condition"] = data["Experimental_Condition"].astype("string").str.strip()
    data["Experiment_Date"] = data["Experiment_Date"].astype("string").str.strip()
    data["Cell death score"] = pd.to_numeric(data["Cell death score"], errors="coerce")
    if data.isna().any(axis=None):
        rows = (data.index[data.isna().any(axis=1)] + 2).tolist()
        fail(f"{path.name} has missing or unreadable values in row(s): {rows}")
    if data.empty:
        fail(f"{path.name} holds no observations")
    bad_dates = sorted({d for d in data["Experiment_Date"]
                        if not pd.Series([d]).str.fullmatch(r"\d{4}-\d{2}-\d{2}").iloc[0]})
    if bad_dates:
        fail(f"dates must be written as YYYY-MM-DD; found: {bad_dates}")
    scores = data["Cell death score"]
    if scores.lt(SCORE_MIN).any() or scores.gt(SCORE_MAX).any():
        fail(f"scores must lie on the {SCORE_MIN} to {SCORE_MAX} scale; "
             f"found {scores.min():g} to {scores.max():g}")
    present = data["Experimental_Condition"].drop_duplicates().tolist()
    conditions = tuple(requested) if requested else tuple(present)
    if len(conditions) != 2:
        fail(f"exactly two conditions are needed; found {len(present)}: {', '.join(present)}")
    unknown = [c for c in conditions if c not in present]
    if unknown:
        fail(f"condition(s) not in {path.name}: {', '.join(unknown)}")
    extra = [c for c in present if c not in conditions]
    if extra:
        fail(f"{path.name} also holds condition(s) that were not selected: {', '.join(extra)}")
    order = data["NLR"].drop_duplicates().tolist()
    data["NLR"] = pd.Categorical(data["NLR"], categories=order, ordered=True)
    return data, conditions


def test_each_receptor(data, conditions):
    rows = []
    for nlr, group in data.groupby("NLR", observed=True, sort=False):
        condition = group["Experimental_Condition"].astype("string")
        first = group.loc[condition == conditions[0], "Cell death score"].to_numpy(float)
        second = group.loc[condition == conditions[1], "Cell death score"].to_numpy(float)
        if len(first) == 0 or len(second) == 0:
            fail(f"{nlr} is missing one of the two conditions")
        result = mannwhitneyu(first, second, alternative="two-sided",
                              method="asymptotic", use_continuity=True)
        pairs = len(first) * len(second)
        rows.append({
            "NLR": str(nlr),
            "Test": TEST_NAME,
            "Date_Stratified": False,
            "n_Replicates": int(group["Experiment_Date"].nunique()),
            "Experimental_Condition_1": conditions[0],
            "n_1": len(first),
            "median_1": float(np.median(first)),
            "Experimental_Condition_2": conditions[1],
            "n_2": len(second),
            "median_2": float(np.median(second)),
            "U_1": float(result.statistic),
            "n_Pairs": pairs,
            "Probability_Superiority_1_over_2": float(result.statistic) / pairs,
            "p_value": float(result.pvalue),
            "significance": significance(float(result.pvalue)),
        })
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    data, conditions = load(args.input_tsv, args.conditions)
    results = test_each_receptor(data, conditions)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Step 8.1 test_cell_death_scores completed: "
          f"observations={len(data)}, receptors={len(results)}, "
          f"conditions={conditions[0]} vs {conditions[1]}, "
          f"dates={data['Experiment_Date'].nunique()}, "
          f"significant_at_0.05={int((results['p_value'] < 0.05).sum())}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
