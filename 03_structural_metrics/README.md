# Structural metrics of the predicted oligomers

Each of the 296 candidate NLR-IDs was modeled three times as a hexamer
(Step 2.5), and the three predictions per receptor are summarized here as six
metrics. Five describe the geometry of the predicted oligomer and one is
AlphaFold's interface confidence. The helper NLR NbNRC2a, whose hexameric
resistosome has been solved by cryo-electron microscopy, was modeled and
scored in the same way as a benchmark, giving 297 scored entries.

| Metric | Meaning | Favorable |
| --- | --- | --- |
| `ipTM_LCB` | Lower confidence bound of the interface predicted template modeling score reported by AlphaFold 3 (referred to as ipTM in the manuscript) | higher |
| `D_APEX` | Mean pairwise distance between the N-terminal methionines at the predicted resistosome apex | lower |
| `SD_THETA_ROT` | Standard deviation of the angles between adjacent protomer centroids around the ring center (`σ_θ ROT`) | lower |
| `S_PROTO` | Radial dispersion of protomer centroids (root-mean-square deviation of their distances from the ring center) | lower |
| `Sum_CONTACTS` | Number of interprotomer contacts (`Σ_CONTACTS`) | higher |
| `BSA_INT_PROTO` | Interprotomer buried surface area (`BSA_INTER-PROTO`) | higher |

Each metric is reported as a confidence bound over the three predictions,
calculated from the replicate mean and sample standard deviation as
`mean ± 1.96 × SD / sqrt(3)`: a lower bound where higher values are favorable,
and an upper bound where lower values are favorable.

The metrics were computed with the pipeline of
[`amiralito/SolNRCH_foldome`](https://github.com/amiralito/SolNRCH_foldome), and
the values below are its output. The metrics are defined in Materials and Methods
of the manuscript.

That pipeline accompanies:

> Toghani, A., Seager, B.A., Sugihara, Y., Roijen, L.-M., Azcue, J.M., Garro, M.,
> Sargolzaei, M., Morianou, I., Harant, A., Gallop, S., Kourelis, J., MacLean, D.,
> Contreras, M.P., Kamoun, S., and Lüdke, D. AI-guided discovery of atypical
> protein assemblies. *bioRxiv* (2026).
> [doi:10.64898/2026.05.03.722499](https://doi.org/10.64898/2026.05.03.722499)

## Step 3.1: Score the predicted oligomers

### Purpose

Score the three predictions of each modeled region with the six metrics above.

### Script and command

Scoring was run with the `SolNRCH_foldome` pipeline cited above, which is
maintained in its own repository and is not duplicated here. Its input is the
AlphaFold 3 output from Step 2.5.

### Outputs

- [`3.1_structural_scores.tsv`](data/results/3.1_structural_scores.tsv): the six
  metrics for the 296 candidates and the NbNRC2a benchmark; 297 rows. These are
  the values published as Table S1 and serve as the input to the ranking.

`ID` identifies a modeled region as `<assembly accession>_6x<protein
accession>_<start>_<end>`, in lowercase with dots replaced by underscores. It
corresponds to the identifiers written in Step 2.1.
