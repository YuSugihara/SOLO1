# Potato virus X infection assay

This workflow contains the fluorescence measurements from the potato virus X
(PVX) infection assay and the statistical tests comparing treatments within each
leaf. It forms the basis of Fig. 4C–D and Tables S6–S7 of the manuscript.

Two reporter viruses were assayed in parallel: PVX expressing GFP (`PVX-GFP`)
and PVX expressing mCherry (`PVX-mCherry`). Each channel is stored and tested
separately.

## Inputs

- [`9.0_pvx_fluorescence.tsv`](data/inputs/9.0_pvx_fluorescence.tsv): 168
  measurements across seven treatments in each of two channels, collected over
  three experiment dates.

The file is in long format, one infiltration site per row, with the columns
`Channel`, `Receptor`, `Leaf`, `Experiment_Date`, `Replicate`, and
`Mean intensity per area`. Dates are written as `YYYY-MM-DD`, and `Replicate`
numbers the three experiments in chronological order. The seven treatments are
the buffer-only control, the empty vector (EV), Rx, the
Pikm-1<sup>Enh</sup>/Pikm-2 pair, and SOLO1<sup>WT</sup>, SOLO1<sup>HIPP43</sup>, and
SOLO1<sup>Enh</sup>. Every treatment was infiltrated into a separate site on
each of the four leaves per experiment (12 leaves per channel), and fluorescence
was recorded three days post inoculation (3 dpi) as the mean intensity per
square centimeter of the infiltrated area.

## Step 9.1: Compare the treatments within each leaf

### Purpose

Test each channel separately. Because every treatment was infiltrated into a
separate site on the same leaves, each of the 21 pairwise comparisons among the
seven treatments is tested within leaves using a two-sided Wilcoxon signed-rank
test on the unnormalized fluorescence values without transformation.

The 21 comparisons within each channel are adjusted together using the Holm
method, and a compact letter display is derived from the adjusted P values at a
family-wise significance level of 0.05: two treatments share a letter when their
adjusted P value is 0.05 or higher.

The output reports the median of each treatment, the median of the paired
differences, and the number of leaves on which the second treatment of a
comparison showed the higher value, so that the direction and magnitude of each
difference can be read alongside the P value.

### Script and command

```bash
python3 09_pvx_infection_assay/scripts/9.1_test_fluorescence_intensity.py \
  --input-tsv 09_pvx_infection_assay/data/inputs/9.0_pvx_fluorescence.tsv \
  --output-summary-tsv 09_pvx_infection_assay/data/results/9.1_pvx_fluorescence_summary.tsv \
  --output-pairwise-tsv 09_pvx_infection_assay/data/results/9.1_pvx_fluorescence_pairwise.tsv
```

`--alpha` sets the family-wise significance level used for the compact letter
display; it defaults to 0.05.

### Outputs

- [`9.1_pvx_fluorescence_summary.tsv`](data/results/9.1_pvx_fluorescence_summary.tsv):
  one row per treatment and channel (14 rows), reporting descriptive statistics
  and the compact letter display.
- [`9.1_pvx_fluorescence_pairwise.tsv`](data/results/9.1_pvx_fluorescence_pairwise.tsv):
  one row per pairwise comparison (42 rows), reporting unadjusted and adjusted P
  values.
