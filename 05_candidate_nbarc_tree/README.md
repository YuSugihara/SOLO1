# Phylogeny of the modelled candidates

This workflow builds the NB-ARC phylogeny of the 296 candidates that were
modelled with AlphaFold 3, and writes the annotation tracks that are displayed
on it. The tree with its tracks is the basis of Fig. 2A.

## Inputs

- [`2.3_af3_targets.fasta`](../02_af3_modeling/data/results/2.3_af3_targets.fasta):
  the 296 modelled candidates, which define the taxon set.
- [`1.2_nlr_id_filtered_nbarc.fasta`](../01_nlr_id_dataset/data/results/1.2_nlr_id_filtered_nbarc.fasta):
  NB-ARC sequences of the discovery set.
- [`refseq_assemblies.tsv`](../01_nlr_id_dataset/data/inputs/refseq_assemblies.tsv):
  species names for the tip labels.
- [`4.1_ranked_scores.tsv`](../04_ranking/data/results/4.1_ranked_scores.tsv):
  the scores, clade memberships and screen outcomes shown as tracks.

## Step 5.1: Collect the NB-ARC sequences

### Purpose

Take the NB-ARC sequence of each modelled candidate and relabel it as
`<species>|<assembly>|<protein>`.

### Script and command

```bash
python3 05_candidate_nbarc_tree/scripts/5.1_prepare_nbarc_alignment_input.py \
  --target-fasta 02_af3_modeling/data/results/2.3_af3_targets.fasta \
  --nbarc-fasta 01_nlr_id_dataset/data/results/1.2_nlr_id_filtered_nbarc.fasta \
  --assemblies-tsv 01_nlr_id_dataset/data/inputs/refseq_assemblies.tsv \
  --output-fasta 05_candidate_nbarc_tree/data/results/5.1_candidate_nbarc.fasta
```

### Outputs

- [`5.1_candidate_nbarc.fasta`](data/results/5.1_candidate_nbarc.fasta): 296
  NB-ARC sequences.

## Step 5.2: Align the NB-ARC sequences

Align with MAFFT (`--localpair --maxiterate 1000`) and write the alignment on
one line per sequence.

```bash
mafft --maxiterate 1000 --localpair \
  05_candidate_nbarc_tree/data/results/5.1_candidate_nbarc.fasta \
  | seqtk seq -l0 \
  > 05_candidate_nbarc_tree/data/results/5.2_candidate_nbarc.mafft.fasta
```

Output: [`5.2_candidate_nbarc.mafft.fasta`](data/results/5.2_candidate_nbarc.mafft.fasta).

## Step 5.3: Trim the alignment

Remove columns that are gaps in more than 90% of the sequences.

```bash
clipkit 05_candidate_nbarc_tree/data/results/5.2_candidate_nbarc.mafft.fasta \
  -m gappy -g 0.9 \
  -o 05_candidate_nbarc_tree/data/results/5.3_candidate_nbarc.mafft.clipkit.fasta
```

Output: [`5.3_candidate_nbarc.mafft.clipkit.fasta`](data/results/5.3_candidate_nbarc.mafft.clipkit.fasta).

## Step 5.4: Infer the maximum-likelihood tree

Let ModelFinder choose the substitution model and support the branches with
1,000 ultrafast bootstrap replicates.

```bash
iqtree -s 05_candidate_nbarc_tree/data/results/5.3_candidate_nbarc.mafft.clipkit.fasta \
  -pre 05_candidate_nbarc_tree/data/results/5.4_candidate_nbarc \
  -st AA -m MFP -bb 1000 -nt 8
```

The model selected by BIC is `Q.mammal+F+R7`.

### Outputs

- [`5.4_candidate_nbarc.treefile`](data/results/5.4_candidate_nbarc.treefile):
  the maximum-likelihood tree with ultrafast bootstrap support.
- [`5.4_candidate_nbarc.contree`](data/results/5.4_candidate_nbarc.contree):
  the bootstrap consensus tree.
- [`5.4_candidate_nbarc.iqtree`](data/results/5.4_candidate_nbarc.iqtree): the
  run report, including model selection.

## Step 5.5: Root the tree at its midpoint

### Purpose

Midpoint-root the maximum-likelihood tree with Biopython.

### Script and command

```bash
python3 05_candidate_nbarc_tree/scripts/5.5_midpoint_root_tree.py \
  --input-tree 05_candidate_nbarc_tree/data/results/5.4_candidate_nbarc.treefile \
  --output-tree 05_candidate_nbarc_tree/data/results/5.5_candidate_nbarc.midpoint_rooted.treefile
```

### Outputs

- [`5.5_candidate_nbarc.midpoint_rooted.treefile`](data/results/5.5_candidate_nbarc.midpoint_rooted.treefile):
  the tree displayed in Fig. 2A.

## Step 5.6: Write the annotation tracks

### Purpose

Write the four iTOL tracks, keyed on the tip labels of the tree: a colour strip
for the sensor clades, a colour strip for the outcome of the functional screen,
a bar track of the apex distance, and a two-field bar track of the interface
confidence split at 0.5. The first two and the interface confidence are
displayed on Fig. 2A; the apex distance is displayed on Fig. S2.

### Script and command

```bash
python3 05_candidate_nbarc_tree/scripts/5.6_make_itol_annotations.py \
  --tree 05_candidate_nbarc_tree/data/results/5.5_candidate_nbarc.midpoint_rooted.treefile \
  --table 04_ranking/data/results/4.1_ranked_scores.tsv \
  --output-dir 05_candidate_nbarc_tree/data/results
```

### Outputs

- [`5.6_itol_sensor_clades.txt`](data/results/5.6_itol_sensor_clades.txt): 165
  tips, 109 Pia-2/MIC1, 42 NRC-S and 14 Pii-2.
- [`5.6_itol_screen_outcome.txt`](data/results/5.6_itol_screen_outcome.txt): the
  eight tested receptors.
- [`5.6_itol_apex_distance.txt`](data/results/5.6_itol_apex_distance.txt) and
  [`5.6_itol_interface_score.txt`](data/results/5.6_itol_interface_score.txt):
  296 tips each.

## Step 5.7: Display the tree

Upload [`5.5_candidate_nbarc.midpoint_rooted.treefile`](data/results/5.5_candidate_nbarc.midpoint_rooted.treefile)
to [iTOL](https://itol.embl.de), then the four annotation files from Step 5.6.
The benchmark NbNRC2a is not part of this tree, which holds the 296 candidates
only, so it carries no track.
