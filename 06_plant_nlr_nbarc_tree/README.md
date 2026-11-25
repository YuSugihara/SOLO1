# Phylogeny of plant NLRs

This workflow builds the NB-ARC phylogeny placing the SOLO1 clade among monocot
and reference plant NLRs (Fig. S4). It combines the NLRs annotated in 22
monocot genomes with the curated reference set RefPlantNLR
([Kourelis et al., 2021](https://doi.org/10.1371/journal.pbio.3001124)) and
colors the tips by the N-terminal domain of each receptor.

## Inputs

- [`6.0_genome_manifest.tsv`](data/inputs/6.0_genome_manifest.tsv): the 22
  assemblies—17 Poaceae genomes as the ingroup and 5 non-Poaceae monocot
  genomes as outgroups—with the taxon prefix used in the sequence identifiers.
- [`6.0_nlrtracker/`](data/inputs/6.0_nlrtracker): the NLRtracker output for
  each assembly, from the dataset deposited by Toghani and Kamoun (2024),
  [doi:10.5281/zenodo.13627395](https://doi.org/10.5281/zenodo.13627395).
- [`6.0_refplantnlr/`](data/inputs/6.0_refplantnlr): the RefPlantNLR sequences
  with their NLRtracker classification. Step 6.2 reads
  `RefPlantNLR_NBARC.fasta` (488 NB-ARC records) and
  `RefPlantNLR_NLRtracker.tsv`. `RefPlantNLR_NLR.fasta` contains the 477
  full-length proteins those NB-ARC records were extracted from; it is read by
  Step 7.14 in [`07_solo1_homologues`](../07_solo1_homologues) to extract
  `NbNRC4c`.

## Step 6.1: Collect the NLRs of the selected genomes

### Purpose

For each assembly, retain proteins that NLRtracker classifies with
`Status = NLR`; degenerate NLRs and every other status are excluded. Write the
full-length and NB-ARC sequences with identifiers that carry the taxon prefix
and end in the simplified domain architecture so that the architecture can be
read back from each tree tip. A protein with multiple NB-ARC regions gets one
record per region, numbered `_NBARC1`, `_NBARC2`, and so on.

The full-length sequences are written alongside the NB-ARC sequences in
`data/results/` and serve as the input to Step 7.1 in
[`07_solo1_homologues`](../07_solo1_homologues).

### Script and command

```bash
python3 06_plant_nlr_nbarc_tree/scripts/6.1_build_nbarc_dataset.py \
  --genome-manifest 06_plant_nlr_nbarc_tree/data/inputs/6.0_genome_manifest.tsv \
  --nlrtracker-dir 06_plant_nlr_nbarc_tree/data/inputs/6.0_nlrtracker \
  --output-full-length-fasta 06_plant_nlr_nbarc_tree/data/results/6.1_selected_genomes_full_length_NLR.fasta \
  --output-nbarc-fasta 06_plant_nlr_nbarc_tree/data/results/6.1_selected_genomes_nbarc.fasta \
  --output-summary-tsv 06_plant_nlr_nbarc_tree/data/results/6.1_selected_genomes_summary.tsv
```

### Outputs

- [`6.1_selected_genomes_nbarc.fasta`](data/results/6.1_selected_genomes_nbarc.fasta):
  16,723 NB-ARC records from 15,828 NLRs; 1,033 records of other status were
  excluded.
- [`6.1_selected_genomes_full_length_NLR.fasta`](data/results/6.1_selected_genomes_full_length_NLR.fasta):
  the same 15,828 full-length NLRs (used by Step 7.1 to extract the SOLO1
  clade).
- [`6.1_selected_genomes_summary.tsv`](data/results/6.1_selected_genomes_summary.tsv):
  per-assembly counts.

## Step 6.2: Annotate the reference NB-ARC identifiers

### Purpose

Append the simplified domain architecture to each RefPlantNLR identifier,
following the same convention as Step 6.1, so that reference and genome
sequences can be colored by the same rule.

### Script and command

```bash
python3 06_plant_nlr_nbarc_tree/scripts/6.2_annotate_refplantnlr_nbarc.py \
  --input-nbarc-fasta 06_plant_nlr_nbarc_tree/data/inputs/6.0_refplantnlr/RefPlantNLR_NBARC.fasta \
  --input-nlrtracker-tsv 06_plant_nlr_nbarc_tree/data/inputs/6.0_refplantnlr/RefPlantNLR_NLRtracker.tsv \
  --output-fasta 06_plant_nlr_nbarc_tree/data/results/6.2_refplantnlr_nbarc.fasta
```

### Outputs

- [`6.2_refplantnlr_nbarc.fasta`](data/results/6.2_refplantnlr_nbarc.fasta): 487
  records; one record of other status was excluded.

## Step 6.3: Merge the two sets and filter by length (250–400 aa)

### Purpose

Combine the two NB-ARC sets, strip parentheses from identifiers (because several
downstream tools interpret them as Newick syntax), and retain sequences between
250 and 400 amino acids in length.

### Script and command

```bash
python3 06_plant_nlr_nbarc_tree/scripts/6.3_merge_and_filter_nbarc.py \
  --input-fasta 06_plant_nlr_nbarc_tree/data/results/6.1_selected_genomes_nbarc.fasta \
                06_plant_nlr_nbarc_tree/data/results/6.2_refplantnlr_nbarc.fasta \
  --output-fasta 06_plant_nlr_nbarc_tree/data/results/6.3_merged_nbarc.len250-400.fasta \
  --report-tsv 06_plant_nlr_nbarc_tree/data/results/6.3_merged_nbarc_length_report.tsv
```

### Outputs

- [`6.3_merged_nbarc.len250-400.fasta`](data/results/6.3_merged_nbarc.len250-400.fasta):
  14,364 sequences of the 17,210 merged; 2,514 were too short and 332 too long.
- [`6.3_merged_nbarc_length_report.tsv`](data/results/6.3_merged_nbarc_length_report.tsv):
  every merged sequence with its length and filter decision.

## Step 6.4: Align the sequences

Align the 14,364 sequences with FAMSA.

```bash
mkdir -p 06_plant_nlr_nbarc_tree/work
famsa -t 14 -refine_mode on \
  06_plant_nlr_nbarc_tree/data/results/6.3_merged_nbarc.len250-400.fasta \
  06_plant_nlr_nbarc_tree/work/6.4_merged_nbarc.famsa.fasta
```

The alignment is written to the working directory and is not distributed in the
repository because it is 70 MB and is regenerated by this single command.

## Step 6.5: Trim the alignment

```bash
clipkit 06_plant_nlr_nbarc_tree/work/6.4_merged_nbarc.famsa.fasta \
  -m gappy -g 0.9 \
  -o 06_plant_nlr_nbarc_tree/data/results/6.5_merged_nbarc.famsa.clipkit.fasta
```

Output: [`6.5_merged_nbarc.famsa.clipkit.fasta`](data/results/6.5_merged_nbarc.famsa.clipkit.fasta).

## Step 6.6: Infer the tree

Infer the tree with VeryFastTree.

```bash
VeryFastTree -threads 14 \
  06_plant_nlr_nbarc_tree/data/results/6.5_merged_nbarc.famsa.clipkit.fasta \
  > 06_plant_nlr_nbarc_tree/data/results/6.6_merged_nbarc.vft.tree
```

Output: [`6.6_merged_nbarc.vft.tree`](data/results/6.6_merged_nbarc.vft.tree), 14,364 tips.

## Step 6.7: Color the tips by N-terminal domain

### Purpose

Read the domain-architecture label from the last underscore-delimited field of
each tip label and color the tip when it contains exactly one of the CC (`C`),
TIR (`T`), or RPW8 (`R`) codes. A CC architecture that also carries one or more
integrated domains (`O`) receives its own color. Tips whose architecture is
absent or ambiguous remain uncolored, and the classification report records why.

### Script and command

```bash
python3 06_plant_nlr_nbarc_tree/scripts/6.7_make_itol_tip_symbols.py \
  --input-tree 06_plant_nlr_nbarc_tree/data/results/6.6_merged_nbarc.vft.tree \
  --output-tip-symbols 06_plant_nlr_nbarc_tree/data/results/6.7_itol_n_terminal_domain.txt \
  --output-report-tsv 06_plant_nlr_nbarc_tree/data/results/6.7_n_terminal_domain_classification.tsv
```

### Outputs

- [`6.7_itol_n_terminal_domain.txt`](data/results/6.7_itol_n_terminal_domain.txt):
  11,522 colored tips, comprising 10,265 CC, 1,071 CC with integrated
  domain(s), 146 TIR, and 40 RPW8.
- [`6.7_n_terminal_domain_classification.tsv`](data/results/6.7_n_terminal_domain_classification.tsv):
  all 14,364 tips with the detected architecture and, for the 2,842 uncolored
  tips, the reason.

## Step 6.8: Label the reference receptors

- [`6.8_itol_reference_labels.txt`](data/results/6.8_itol_reference_labels.txt):
  external iTOL labels in bold black for the characterized reference receptors
  shown in the figure.

Steps 6.7–6.8 provide optional iTOL annotation files; Fig. S4 itself is plotted
directly in R in Step 6.9.

## Step 6.9: Plot the plant NLR phylogeny (Fig. S4)

### Purpose

Root the 14,364-tip NB-ARC phylogeny on the RNL/TNL outgroup clade and render it
as an open-fan tree directly in R using `ggtree`, with two concentric outer
tracks: (i) N-terminal domain classes from Step 6.7 (`CC with ID`,
`CC without ID`, `TIR`, `RPW8`, `Unclassified`) and (ii) focal clade and ID
annotations. These annotations cover the monophyletic `SOLO1 clade`
(`n = 136`, SH-aLRT = `1.000`), `DDE` (`176`), and `BED` (`160`), together with
`MIC1` (`n = 582`, SH-aLRT = `1.000`, `66.0%` CC with ID), `MIC2` (`n = 367`,
`169/176` DDE), `MIC3` (`n = 232`, SH-aLRT = `0.987`, `160/160` BED), and 21
reference NLRs connected to their tree tips by dotted radial guides.

### Script and command

```bash
Rscript 06_plant_nlr_nbarc_tree/scripts/6.9_plot_plant_nlr_phylogeny.R \
  --tree 06_plant_nlr_nbarc_tree/data/results/6.6_merged_nbarc.vft.tree \
  --classification 06_plant_nlr_nbarc_tree/data/results/6.7_n_terminal_domain_classification.tsv \
  --solo1-ids 07_solo1_homologues/data/inputs/7.0_solo1_clade_ids.txt \
  --id-gff 01_nlr_id_dataset/data/results/1.1_nlr_id_candidates.gff \
  --nlrtracker-table 01_nlr_id_dataset/data/inputs/refseq_nlrtracker_tables.tsv.gz \
  --output-prefix 06_plant_nlr_nbarc_tree/data/results/6.9_plant_nlr_phylogeny
```

### Outputs

- [`6.9_plant_nlr_phylogeny.pdf`](data/results/6.9_plant_nlr_phylogeny.pdf):
  vector figure (Type 1 `Helvetica`) for Fig. S4.
- [`6.9_plant_nlr_phylogeny.png`](data/results/6.9_plant_nlr_phylogeny.png):
  300 DPI preview PNG.
