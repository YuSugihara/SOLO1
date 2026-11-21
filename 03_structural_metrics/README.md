# Structural metrics of the predicted oligomers

Each of the 296 candidate NLR-IDs was modelled three times as a hexamer
(Step 2.4), and the three predictions of a receptor are summarised here as six
metrics. Five describe the geometry of the predicted oligomer and one is
AlphaFold's own interface confidence. The helper NLR NbNRC2a, whose hexameric
resistosome has been solved by cryo-electron microscopy, was modelled and
scored in the same way and serves as the benchmark, giving 297 scored entries.

| Metric | Meaning | Favourable |
| --- | --- | --- |
| `ipTM_LCB` | Interface predicted template modelling score reported by AlphaFold 3 | higher |
| `D_APEX` | Average distance between the N-terminal methionines at the predicted resistosome apex | lower |
| `SD_THETA_ROT` | Dispersion of the rotational angle between protomers | lower |
| `S_PROTO` | Radial dispersion of protomer centroids, the root-mean-square deviation of their distances from the ring centre | lower |
| `Sum_CONTACTS` | Number of contacts between protomers | higher |
| `BSA_INT_PROTO` | Surface area buried between protomers | higher |

Each metric is reported as a confidence bound over the three predictions,
calculated from the replicate mean and sample standard deviation as
`mean ± 1.96 × SD / sqrt(3)`: a lower bound where higher values are favourable,
and an upper bound where lower values are favourable.

The calculations are described in the Methods of the paper.

## Outputs

- [`3.1_structural_scores.tsv`](data/results/3.1_structural_scores.tsv): the six
  metrics for the 296 candidates and the NbNRC2a benchmark; 297 rows. These are
  the values published as Table S1 and are the input to the ranking.

`ID` identifies a modelled region as `<assembly accession>_6x<protein
accession>_<start>_<end>`, in lower case with dots replaced by underscores. It
corresponds to the identifiers written in Step 2.1.
