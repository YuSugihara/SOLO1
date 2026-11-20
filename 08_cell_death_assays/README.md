# Cell death assays

This workflow holds the cell death scores of the transient expression assays and
the test that compares the two conditions of each assay. It is the basis of
Fig. 3B and fig. S7, and of tables S2 to S5 of the paper.

The two assays are held separately. The first compares each wild-type receptor
with its MHD-motif aspartate-to-valine (DV) variant. The second compares an
engineered receptor co-expressed with its matched ligand against the same
receptor co-expressed with the mismatched ligand.

## Inputs

- [`8.0_wt_vs_dv_scores.tsv`](data/inputs/8.0_wt_vs_dv_scores.tsv): 180 scores,
  six SOLO1 homologues, each as wild type (WT) and as its DV variant, over five
  experiment dates.
- [`8.0_gfp_vs_pwl2_scores.tsv`](data/inputs/8.0_gfp_vs_pwl2_scores.tsv): 104
  scores, four engineered receptors, each with GFP and with PWL2, over three
  experiment dates.

Both files are in long format, one observation per row, with the columns
`Comparison`, `NLR`, `Experimental_Condition`, `Experiment_Date` and
`Cell death score`. Dates are written as `YYYY-MM-DD`. Cell death was scored
visually on a 0 to 7 scale, where 0 is no visible symptom and 7 is confluent
necrosis across the infiltrated area.

## Step 8.1: Compare the two conditions of each assay

### Purpose

Test each receptor on its own. The two conditions are compared with a two-sided
Mann-Whitney U test on the pooled observations, using the asymptotic
approximation with the tie correction and the continuity correction. The scores
are an ordinal scale with many ties, so the exact test does not apply.

The experiment date is read, checked and counted, but it takes no part in the
test: the observations of all dates are pooled.

The output carries the sample size and the median of each condition, the U
statistic of the first condition and the probability of superiority alongside
the P value.

### Script and command

```bash
python3 08_cell_death_assays/scripts/8.1_test_cell_death_scores.py \
  --input-tsv 08_cell_death_assays/data/inputs/8.0_wt_vs_dv_scores.tsv \
  --output-tsv 08_cell_death_assays/data/results/8.1_wt_vs_dv_mann_whitney.tsv

python3 08_cell_death_assays/scripts/8.1_test_cell_death_scores.py \
  --input-tsv 08_cell_death_assays/data/inputs/8.0_gfp_vs_pwl2_scores.tsv \
  --output-tsv 08_cell_death_assays/data/results/8.1_gfp_vs_pwl2_mann_whitney.tsv
```

### Outputs

- [`8.1_wt_vs_dv_mann_whitney.tsv`](data/results/8.1_wt_vs_dv_mann_whitney.tsv):
  one row per receptor, six rows.
- [`8.1_gfp_vs_pwl2_mann_whitney.tsv`](data/results/8.1_gfp_vs_pwl2_mann_whitney.tsv):
  one row per receptor, four rows.

## Step 8.2: Draw the cell death score panel

### Purpose

Draw Fig. 3B and fig. S7. Each receptor occupies one block on a shared y axis
of scores. Within a block the two conditions are drawn side by side, and within
a condition each experiment date is drawn as its own column. Observations that
share a score within a column are collapsed into one point whose area is
proportional to the number of observations. The P value of each receptor is
read from the table written by Step 8.1.

### Script and command

```bash
python3 08_cell_death_assays/scripts/8.2_plot_cell_death_scores.py \
  --input-tsv 08_cell_death_assays/data/inputs/8.0_wt_vs_dv_scores.tsv \
  --statistics-tsv 08_cell_death_assays/data/results/8.1_wt_vs_dv_mann_whitney.tsv \
  --output-prefix 08_cell_death_assays/data/results/8.2_wt_vs_dv_scores

python3 08_cell_death_assays/scripts/8.2_plot_cell_death_scores.py \
  --input-tsv 08_cell_death_assays/data/inputs/8.0_gfp_vs_pwl2_scores.tsv \
  --statistics-tsv 08_cell_death_assays/data/results/8.1_gfp_vs_pwl2_mann_whitney.tsv \
  --output-prefix 08_cell_death_assays/data/results/8.2_gfp_vs_pwl2_scores \
  --colours "#2C8C8C" "#B5367A"
```

### Outputs

- [`8.2_wt_vs_dv_scores.pdf`](data/results/8.2_wt_vs_dv_scores.pdf),
  [`.svg`](data/results/8.2_wt_vs_dv_scores.svg) and
  [`.png`](data/results/8.2_wt_vs_dv_scores.png): the panel published as
  fig. S7.
- [`8.2_gfp_vs_pwl2_scores.pdf`](data/results/8.2_gfp_vs_pwl2_scores.pdf),
  [`.svg`](data/results/8.2_gfp_vs_pwl2_scores.svg) and
  [`.png`](data/results/8.2_gfp_vs_pwl2_scores.png): the panel published as
  Fig. 3B.

Text is kept as text in the vector formats.
