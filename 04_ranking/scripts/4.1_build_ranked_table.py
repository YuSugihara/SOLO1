#!/usr/bin/env python3
"""Step 4.1 script: rank the candidates on their structural metrics.

Each of the six metrics is turned into a favourable rank across all scored
entries, so that rank 1 is the most favourable value whichever direction of the
metric is better. The entries are then ordered on a weighted combination of two
of those ranks, the interface confidence and the apex distance, which is the
ranking reported in the paper. The other four ranks are carried through for
comparison and take no part in the ordering.

Entries are also annotated with the sensor clade they belong to, and with the
outcome of the functional screen for the receptors that were tested.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCORES = ["ipTM_LCB", "BSA_INT_PROTO", "Sum_CONTACTS",
          "S_PROTO", "SD_THETA_ROT", "D_APEX"]
LOWER_IS_BETTER = {"S_PROTO", "SD_THETA_ROT", "D_APEX"}
CLADE_FILES = {"Pii-2": "pii2_clade.txt", "NRC-S": "nrcs_clade.txt",
               "Pia-2/MIC1": "mic1_clade.txt"}
NOT_APPLICABLE = "N/A"
AUTOACTIVE = {"Autoactive", "Strong autoactive"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scores", required=True, type=Path,
                   help="Structural metrics, one row per scored entry.")
    p.add_argument("--annotations", required=True, type=Path,
                   help="Identifiers, species and functional-screen outcomes.")
    p.add_argument("--clade-dir", required=True, type=Path,
                   help="Directory holding the three sensor-clade lists.")
    p.add_argument("--output-tsv", required=True, type=Path,
                   help="Ranked and annotated table.")
    p.add_argument("--reference-id", default="NbNRC2a",
                   help="Entry used as the hexamer benchmark (default: %(default)s).")
    p.add_argument("--iptm-weight", type=float, default=0.70,
                   help="Weight of the ipTM rank (default: %(default)s).")
    p.add_argument("--apex-weight", type=float, default=0.30,
                   help="Weight of the apex-distance rank (default: %(default)s).")
    return p.parse_args()


def fail(message):
    sys.exit(f"error: {message}")


def load_inputs(scores_path, annotations_path):
    scores = pd.read_csv(scores_path, sep="\t", keep_default_na=False)
    annotations = pd.read_csv(annotations_path, sep="\t", keep_default_na=False)
    for frame, path, needed in ((scores, scores_path, {"ID", *SCORES}),
                                (annotations, annotations_path,
                                 {"ID", "protein_ID", "species", "Cell_death",
                                  "WT_domesticated_autoactivity"})):
        missing = sorted(needed.difference(frame.columns))
        if missing:
            fail(f"{path} is missing columns: {missing}")
        if frame["ID"].duplicated().any():
            fail(f"{path} contains duplicate IDs")
    if set(scores["ID"]) != set(annotations["ID"]):
        fail("the two inputs cover different entries")
    data = annotations.merge(scores, on="ID", validate="one_to_one")
    for score in SCORES:
        data[score] = pd.to_numeric(data[score], errors="raise")
    if not np.isfinite(data[SCORES].to_numpy()).all():
        fail("every structural metric must be a finite number")
    if data["protein_ID"].duplicated().any():
        fail("protein_ID must be unique")
    return data


def classify_screen(data):
    """Derive the screen outcome from the two recorded observations.

    Receptors that were not part of the screen reported here carry no
    observation and stay unannotated.
    """
    outcome = pd.Series(NOT_APPLICABLE, index=data.index, dtype="object")
    autoactive = data["WT_domesticated_autoactivity"].isin(AUTOACTIVE)
    outcome.loc[data["Cell_death"].eq("No")] = "No cell death"
    outcome.loc[data["Cell_death"].eq("Yes") & autoactive] = "Autoactive"
    outcome.loc[data["Cell_death"].eq("Yes") & ~autoactive] = "DV-dependent cell death"
    return outcome


def rank_entries(data, iptm_weight, apex_weight):
    ranked = data.copy()
    total = len(ranked)
    for score in SCORES:
        ranked[f"{score}_favourable_rank_{total}"] = ranked[score].rank(
            method="average", ascending=score in LOWER_IS_BETTER)
    ranked["ipTM_LCB_weight"] = iptm_weight
    ranked["D_APEX_weight"] = apex_weight
    # Rounded before ordering: the ranks are multiples of 0.5, so the exact
    # scores are well separated, and without rounding two equal scores can
    # differ by ~1e-14 of floating-point noise and fail to tie.
    ranked["weighted_priority_score"] = (
        iptm_weight * ranked[f"ipTM_LCB_favourable_rank_{total}"]
        + apex_weight * ranked[f"D_APEX_favourable_rank_{total}"]).round(6)
    # Entries that share a composite score are ordered by ipTM, the more
    # heavily weighted of the two metrics; ID breaks any remaining tie so that
    # the order is reproducible.
    ranked = ranked.sort_values(["weighted_priority_score", "ipTM_LCB", "ID"],
                                ascending=[True, False, True],
                                kind="mergesort").reset_index(drop=True)
    ranked["weighted_rank"] = np.arange(1, total + 1)
    return ranked


def add_sensor_clades(data, clade_dir):
    annotated = data.copy()
    annotated["sensor_clade"] = NOT_APPLICABLE
    available = set(annotated["protein_ID"])
    assigned = set()
    for label, filename in CLADE_FILES.items():
        path = clade_dir / filename
        if not path.exists():
            fail(f"{path} does not exist")
        accessions = {line.strip().split("|")[-1]
                      for line in path.read_text(encoding="utf-8").splitlines()
                      if line.strip()}
        overlap = assigned & accessions
        if overlap:
            fail(f"{filename} shares members with another clade: {sorted(overlap)[:5]}")
        missing = accessions - available
        if missing:
            fail(f"{filename} lists entries absent from the scores: {sorted(missing)[:5]}")
        annotated.loc[annotated["protein_ID"].isin(accessions), "sensor_clade"] = label
        assigned |= accessions
    annotated["is_sensor_clade_member"] = annotated["sensor_clade"].ne(NOT_APPLICABLE)
    return annotated


def validate(data, reference_id, total):
    if len(data) != total or data["protein_ID"].nunique() != total:
        fail(f"the table must hold {total} entries with unique protein IDs")
    if data["weighted_rank"].tolist() != list(range(1, total + 1)):
        fail(f"weighted ranks must run from 1 to {total} without gaps")
    # Every entry gets a rank of its own, so a shared composite score must be
    # settled by ipTM. If two entries also share their ipTM the order would
    # fall through to the identifier, which carries no meaning, so stop here.
    for score, group in data.groupby("weighted_priority_score"):
        if len(group) > 1 and group["ipTM_LCB"].nunique() != len(group):
            fail(f"composite score {score} is shared by entries with the same "
                 f"ipTM, so their order would be arbitrary: "
                 f"{sorted(group['protein_ID'])}")
    references = set(data.loc[data["is_reference"], "protein_ID"])
    if references != {reference_id}:
        fail(f"{reference_id} must be the only benchmark; found {sorted(references)}")


def main():
    args = parse_args()
    data = load_inputs(args.scores, args.annotations)
    total = len(data)

    data["is_reference"] = data["protein_ID"].str.lower().eq(args.reference_id.lower())
    if int(data["is_reference"].sum()) != 1:
        fail(f"expected exactly one entry for {args.reference_id}")
    data["screening_class"] = classify_screen(data)
    data.loc[data["is_reference"], "screening_class"] = NOT_APPLICABLE

    ranked = rank_entries(data, args.iptm_weight, args.apex_weight)
    annotated = add_sensor_clades(ranked, args.clade_dir)
    validate(annotated, args.reference_id, total)

    leading = ["weighted_rank", "weighted_priority_score",
               "ipTM_LCB_weight", "D_APEX_weight",
               *[f"{s}_favourable_rank_{total}" for s in SCORES],
               "is_reference", "sensor_clade", "is_sensor_clade_member",
               "screening_class"]
    rest = [c for c in annotated.columns if c not in leading]
    final = annotated[leading + rest]
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.output_tsv, sep="\t", index=False, float_format="%.10g")

    counts = final.loc[~final["is_reference"], "screening_class"].value_counts().to_dict()
    clades = final["sensor_clade"].value_counts().to_dict()
    reference_rank = int(final.loc[final["is_reference"], "weighted_rank"].iloc[0])
    print(f"Step 4.1 build_ranked_table completed: entries={total}, "
          f"{args.reference_id}_rank={reference_rank}, screen={counts}, clades={clades}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
