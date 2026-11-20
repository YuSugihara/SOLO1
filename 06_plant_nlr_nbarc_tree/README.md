# Phylogeny of plant NLRs

This workflow builds the NB-ARC phylogeny that places the SOLO1 clade among
plant NLRs as a whole, which is the basis of Fig. S4. It combines the NLRs
annotated in 22 monocot genomes with the curated reference set RefPlantNLR
([Kourelis et al., 2021](https://doi.org/10.1371/journal.pbio.3001124)), and
colours the tips by the N-terminal domain of each receptor.

## Inputs

- [`6.0_genome_manifest.tsv`](data/inputs/6.0_genome_manifest.tsv): the 22
  assemblies, 17 Poaceae as the ingroup and 5 further monocots as outgroups,
  with the taxon prefix used in the sequence identifiers.
- [`6.0_nlrtracker/`](data/inputs/6.0_nlrtracker): the NLRtracker output of each
  assembly, from the dataset deposited by Toghani and Kamoun (2024),
  [doi:10.5281/zenodo.13627395](https://doi.org/10.5281/zenodo.13627395).
- [`6.0_refplantnlr/`](data/inputs/6.0_refplantnlr): the RefPlantNLR sequences
  with their NLRtracker classification.

## Step 6.1: Collect the NLRs of the selected genomes

### Purpose

For each assembly, retain the proteins that NLRtracker classifies with
`Status = NLR`; degenerate NLRs and every other status are excluded. Write the
full-length and NB-ARC sequences with identifiers that carry the taxon prefix
and end in the simplified domain architecture, so that the architecture can be
read back from a tree tip. A protein with several NB-ARC regions gets one
record per region, numbered `_NBARC1`, `_NBARC2` and so on.

The full-length sequences are large and are written to a working directory
rather than to the results. They are the input to Step 7.1 in
[`07_solo1_homologues`](../07_solo1_homologues).

### Script and command

```bash
python3 06_plant_nlr_nbarc_tree/scripts/6.1_build_nbarc_dataset.py \
  --genome-manifest 06_plant_nlr_nbarc_tree/data/inputs/6.0_genome_manifest.tsv \
  --nlrtracker-dir 06_plant_nlr_nbarc_tree/data/inputs/6.0_nlrtracker \
  --output-full-length-fasta 06_plant_nlr_nbarc_tree/work/6.1_selected_genomes_full_length_NLR.fasta \
  --output-nbarc-fasta 06_plant_nlr_nbarc_tree/data/results/6.1_selected_genomes_nbarc.fasta \
  --output-summary-tsv 06_plant_nlr_nbarc_tree/data/results/6.1_selected_genomes_summary.tsv
```

### Outputs

- [`6.1_selected_genomes_nbarc.fasta`](data/results/6.1_selected_genomes_nbarc.fasta):
  16,723 NB-ARC records from 15,828 NLRs; 1,033 records of other status were
  excluded.
- [`6.1_selected_genomes_summary.tsv`](data/results/6.1_selected_genomes_summary.tsv):
  per-assembly counts.

## Step 6.2: Annotate the reference NB-ARC identifiers

### Purpose

Append the simplified domain architecture to each RefPlantNLR identifier, in the
same convention as Step 6.1, so that reference and genome sequences can be
coloured by the same rule.

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

## Step 6.3: Merge the two sets and keep the well-sized sequences

### Purpose

Combine the two NB-ARC sets, strip parentheses from identifiers, because
several downstream tools treat them as Newick syntax, and keep sequences
between 250 and 400 residues.

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
  every merged sequence with its length and the decision.

## Step 6.4: Align the sequences

Align the 14,364 sequences with FAMSA.

```bash
famsa -t 14 -refine_mode on \
  06_plant_nlr_nbarc_tree/data/results/6.3_merged_nbarc.len250-400.fasta \
  06_plant_nlr_nbarc_tree/work/6.4_merged_nbarc.famsa.fasta
```

The alignment is written to the working directory and is not distributed: it is
70 MB and this one command regenerates it.

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

## Step 6.7: Colour the tips by N-terminal domain

### Purpose

Read the architecture from the last underscore-delimited field of each tip
label and colour the tip when exactly one of the CC, TIR and RPW8 codes is
present. A CC architecture that also carries an integrated domain gets a colour
of its own. Tips whose architecture is absent or ambiguous are left uncoloured,
and the report records why.

### Script and command

```bash
python3 06_plant_nlr_nbarc_tree/scripts/6.7_make_itol_tip_symbols.py \
  --input-tree 06_plant_nlr_nbarc_tree/data/results/6.6_merged_nbarc.vft.tree \
  --output-tip-symbols 06_plant_nlr_nbarc_tree/data/results/6.7_itol_n_terminal_domain.txt \
  --output-report-tsv 06_plant_nlr_nbarc_tree/data/results/6.7_n_terminal_domain_classification.tsv
```

### Outputs

- [`6.7_itol_n_terminal_domain.txt`](data/results/6.7_itol_n_terminal_domain.txt):
  11,522 coloured tips, of which 10,265 CC, 1,071 CC with an integrated domain,
  146 TIR and 40 RPW8.
- [`6.7_n_terminal_domain_classification.tsv`](data/results/6.7_n_terminal_domain_classification.tsv):
  all 14,364 tips with the architecture detected and, for the 2,842 uncoloured
  tips, the reason.

## Step 6.8: Label the reference receptors

- [`6.8_itol_reference_labels.txt`](data/results/6.8_itol_reference_labels.txt):
  external labels for the characterised receptors named in the figure, drawn in
  bold black.

## Step 6.9: Display the tree

Upload [`6.6_merged_nbarc.vft.tree`](data/results/6.6_merged_nbarc.vft.tree) to
[iTOL](https://itol.embl.de), then the tip symbols of Step 6.7 and the labels of
Step 6.8. Enable dashed lines between leaves and external datasets in the
display settings.
