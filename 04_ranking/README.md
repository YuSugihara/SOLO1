# Structural ranking of the candidates

This workflow turns the six structural metrics into the ranking reported in the
paper, and annotates each entry with its sensor clade and, where the receptor
was tested, the outcome of the functional screen.

## Inputs

- [`3.1_structural_scores.tsv`](../03_structural_metrics/data/results/3.1_structural_scores.tsv):
  the six metrics for the 296 candidates and the NbNRC2a benchmark.
- [`4.0_candidate_annotations.tsv`](data/inputs/4.0_candidate_annotations.tsv):
  for each entry, its assembly and protein accession, species, the modelled
  region, the submitted sequence, and the observations from the functional
  screen. `Cell_death` and `WT_domesticated_autoactivity` are recorded for the
  eight receptors screened in this study and are `N/A` for the rest.
- [`4.0_sensor_clades/`](data/inputs/4.0_sensor_clades): the members of the
  Pii-2, NRC-S and Pia-2/MIC1 sensor clades, 14, 42 and 109 accessions
  respectively.

## Step 4.1: Rank the candidates and annotate them

### Purpose

Rank all 297 entries on each metric so that rank 1 is the most favourable value,
taking the direction of each metric into account: higher is better for
`ipTM_LCB`, `Sum_CONTACTS` and `BSA_INT_PROTO`, lower is better for `D_APEX`,
`SD_THETA_ROT` and `S_PROTO`. Tied values share the average of the ranks they
span.

Order the entries on a weighted combination of two of those ranks,

```text
weighted_priority_score = 0.70 x ipTM_LCB rank + 0.30 x D_APEX rank
```

where a smaller score ranks first. Entries that share a composite score are
ordered by ipTM, the more heavily weighted of the two metrics, with the
identifier breaking any remaining tie so that the order is reproducible. The
other four ranks are reported for comparison and take no part in the ordering.

Derive the screen outcome from the two recorded observations: a receptor whose
wild-type form is already autoactive is `Autoactive`, one that gives cell death
only after the MHD-motif substitution is `DV-dependent cell death`, and one that
gives none is `No cell death`. Receptors that were not part of this screen carry
no observation and stay unannotated.

The script stops with an error unless the ranks run from 1 to 297 without gaps,
every shared composite score is settled by ipTM rather than falling through to
the identifier, the clade lists are disjoint and all their members are present,
and the benchmark is the only reference entry.

### Script and command

```bash
python3 04_ranking/scripts/4.1_build_ranked_table.py \
  --scores 03_structural_metrics/data/results/3.1_structural_scores.tsv \
  --annotations 04_ranking/data/inputs/4.0_candidate_annotations.tsv \
  --clade-dir 04_ranking/data/inputs/4.0_sensor_clades \
  --output-tsv 04_ranking/data/results/4.1_ranked_scores.tsv
```

### Outputs

- [`4.1_ranked_scores.tsv`](data/results/4.1_ranked_scores.tsv): all 297 entries
  with their six favourable ranks, the weighted ranking, the sensor clade and
  the screen outcome, together with the two weights and the reference and
  clade-membership flags. Table S1 is drawn from this table.

## Step 4.2: Draw the structural ranking panel

### Purpose

Draw Fig. 1B. The entries are placed along the x axis in the order of the
weighted ranking, with the apex distance drawn as a bar track above and the
interface confidence below. The members of each sensor clade are marked in their
own lane underneath, the benchmark is marked with a dashed line, and the region
that ranks above it is shaded.

### Script and command

```bash
python3 04_ranking/scripts/4.2_plot_structural_ranking.py \
  --table 04_ranking/data/results/4.1_ranked_scores.tsv \
  --output-prefix 04_ranking/data/results/4.2_structural_ranking
```

### Outputs

- [`4.2_structural_ranking.pdf`](data/results/4.2_structural_ranking.pdf),
  [`.svg`](data/results/4.2_structural_ranking.svg) and
  [`.png`](data/results/4.2_structural_ranking.png): the panel published as
  Fig. 1B. Text is kept as text in the vector formats.

The figure of the paper places this panel below a drawing of the workflow,
which is not generated here.
