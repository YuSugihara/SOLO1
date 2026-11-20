# SOLO1

Single-gene NLR immune receptor scaffolds for made-to-order plant disease
resistance.

This repository holds the code and data needed to reproduce the computational
analyses and figures of the SOLO1 study. The study combines AlphaFold 3
structural screening of plant NLR proteins carrying integrated domains
(NLR-IDs) with functional assays in *Nicotiana benthamiana* to identify SOLO1,
a conserved grass NLR clade whose members trigger immune cell death as single
receptors, and to engineer SOLO1 with heterologous ligand-binding modules.

## Repository layout

Each directory is a self-contained analysis with its own `README.md`,
`scripts/` and `data/`. Steps are numbered `<directory>.<step>`, and a script
and the files it writes share that number. Run every command from the root of
the repository.

| Directory | Analysis | Figures |
| --- | --- | --- |
| [`01_nlr_id_dataset`](01_nlr_id_dataset) | Extract CC-NLRs carrying an integrated domain from 180 annotated plant proteomes | Fig. 1A |
| [`02_af3_modeling`](02_af3_modeling) | Prepare the modelled regions and the AlphaFold 3 jobs | Fig. 1A |
| [`03_structural_metrics`](03_structural_metrics) | The six metrics scored on the predicted oligomers | Table S1 |
| [`04_ranking`](04_ranking) | Rank the candidates and draw the ranking panel | Fig. 1B, Table S1 |
| [`05_candidate_nbarc_tree`](05_candidate_nbarc_tree) | NB-ARC phylogeny of the modelled candidates | Fig. 2A, Fig. S2 |
| [`06_plant_nlr_nbarc_tree`](06_plant_nlr_nbarc_tree) | NB-ARC phylogeny of plant NLRs as a whole | Fig. S4 |
| [`07_solo1_homologues`](07_solo1_homologues) | Phylogeny, domain architecture and alignment of the SOLO1 clade | Fig. 2B, Fig. S5, Fig. S6 |
| [`08_cell_death_assays`](08_cell_death_assays) | Cell death scores of the transient expression assays and the test that compares them | Fig. 3B, Fig. S7, Tables S2 to S5 |

Several figures are assembled in [iTOL](https://itol.embl.de) from a tree and
the annotation tracks written here; each analysis says which files to upload.

## Data sources

### Annotated plant proteomes

The NLR-ID discovery set and the genome-wide phylogeny start from the NLRs
annotated with NLRtracker in 180 RefSeq reference plant proteomes. Those
annotations were produced and deposited by Toghani and Kamoun (2024),
[doi:10.5281/zenodo.13627395](https://doi.org/10.5281/zenodo.13627395); the
proteomes themselves are the RefSeq annotations held at NCBI. The deposited
dataset is not redistributed in full. The classification tables of
all 180 assemblies and the per-assembly files needed by
[`06_plant_nlr_nbarc_tree`](06_plant_nlr_nbarc_tree) are included, so only
[Step 1.1](01_nlr_id_dataset) needs the original archive.

### Reference NLRs

The curated reference set RefPlantNLR comes from Kourelis et al. (2021),
[doi:10.1371/journal.pbio.3001124](https://doi.org/10.1371/journal.pbio.3001124).

### Structure predictions

The AlphaFold 3 predictions are too large to distribute through GitHub and are
deposited separately. The scores derived from them are held in
[`03_structural_metrics`](03_structural_metrics), and the submitted jobs in
[`02_af3_modeling`](02_af3_modeling), so the ranking and every figure that rests
on it can be reproduced from this repository alone.

## Software

| Tool | Version |
| --- | --- |
| [Python](https://www.python.org) | 3.10.12 |
| [NumPy](https://numpy.org) | 2.2.6 |
| [pandas](https://pandas.pydata.org) | 2.3.3 |
| [Matplotlib](https://matplotlib.org) | 3.10.9 |
| [SciPy](https://scipy.org) | 1.15.3 |
| [Biopython](https://biopython.org) | 1.88 |
| [CD-HIT](https://sites.google.com/view/cd-hit) | 4.8.1 |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | 7.526 |
| [FAMSA](https://github.com/refresh-bio/FAMSA) | 2.4.1 |
| [ClipKIT](https://jlsteenwyk.com/ClipKIT/) | 2.4.1 |
| [seqtk](https://github.com/lh3/seqtk) | 1.5-r133 |
| [IQ-TREE](http://www.iqtree.org) | 2.3.6 (Step 5.4), 3.0.1 (Step 7.6) |
| [VeryFastTree](https://github.com/citiususc/veryfasttree) | 4.0.5 |
| [TaxonKit](https://bioinf.shenwei.me/taxonkit/) | 0.20.0 |
| [AlphaFold Server](https://alphafoldserver.com) | AlphaFold 3 |

## License

This repository is released under the
[Creative Commons Attribution 4.0 International licence](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0). See [`LICENSE`](LICENSE).

## Citation

Azcue, J.M., Posbeyikian, A., Toghani, A., Chiang, B.-J., Youles, M.,
Lüdke, D., Xi, Y., Kamoun, S., and Sugihara, Y. Single-gene immune receptors
confer made-to-order plant disease resistance. *Submitted*.
