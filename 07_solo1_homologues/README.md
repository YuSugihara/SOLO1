# SOLO1 homologues

This workflow describes the SOLO1 clade: which receptors belong to it, how they
are related, how their domain architecture varies, and how their sequences
align. It produces the tree and the tracks of Fig. S5A and the alignment of
Fig. S6, and the tree that Fig. 2B is drawn from.

## Inputs

- [`7.0_solo1_clade_ids.txt`](data/inputs/7.0_solo1_clade_ids.txt): the 131
  members of the clade. They were read off the plant-wide NLR tree of
  [`06_plant_nlr_nbarc_tree`](../06_plant_nlr_nbarc_tree) by inspecting the clade in a tree
  viewer, so this list is a curated input rather than a computed one.
- [`7.0_seqid_to_assembly.tsv`](data/inputs/7.0_seqid_to_assembly.tsv): the
  assembly each sequence comes from.
- [`7.0_exon_boundaries.tsv`](data/inputs/7.0_exon_boundaries.tsv): exon
  boundaries of the clade members, taken from their RefSeq GFF3 annotations.
- [`7.0_nlrtracker_domains/`](data/inputs/7.0_nlrtracker_domains): the NLRtracker
  domain table of each clade member, one file per assembly, restricted to these
  proteins.

Step 7.1 draws on the full-length NLR set of
[`06_plant_nlr_nbarc_tree`](../06_plant_nlr_nbarc_tree), so run Step 6.1 first.

## Step 7.1: Pull the clade out of the full-length NLR set

```bash
python3 07_solo1_homologues/scripts/7.1_extract_clade_sequences.py \
  --id-list 07_solo1_homologues/data/inputs/7.0_solo1_clade_ids.txt \
  --input-fasta 06_plant_nlr_nbarc_tree/work/6.1_selected_genomes_full_length_NLR.fasta \
  --output-fasta 07_solo1_homologues/data/results/7.1_solo1_clade_full_length.fasta
```

Output: [`7.1_solo1_clade_full_length.fasta`](data/results/7.1_solo1_clade_full_length.fasta), 131 sequences.

## Step 7.2: Remove identical splicing variants

### Purpose

RefSeq annotates several transcripts for some loci, which can yield identical
proteins. Within an assembly, group records that share an identical amino acid
sequence and whose genomic intervals overlap, and keep the record with the
lowest `XP_` accession number from each group.

### Script and command

```bash
python3 07_solo1_homologues/scripts/7.2_remove_identical_variants.py \
  --input-fasta 07_solo1_homologues/data/results/7.1_solo1_clade_full_length.fasta \
  --seqid-assembly-tsv 07_solo1_homologues/data/inputs/7.0_seqid_to_assembly.tsv \
  --exon-boundaries-tsv 07_solo1_homologues/data/inputs/7.0_exon_boundaries.tsv \
  --output-fasta 07_solo1_homologues/data/results/7.2_solo1_clade_nonredundant.fasta \
  --output-removed-tsv 07_solo1_homologues/data/results/7.2_removed_identical_variants.tsv \
  --output-kept-tsv 07_solo1_homologues/data/results/7.2_kept_records.tsv
```

### Outputs

- [`7.2_solo1_clade_nonredundant.fasta`](data/results/7.2_solo1_clade_nonredundant.fasta): 96 sequences; 35 records were removed.
- [`7.2_removed_identical_variants.tsv`](data/results/7.2_removed_identical_variants.tsv) and
  [`7.2_kept_records.tsv`](data/results/7.2_kept_records.tsv): what was removed and what was kept.

## Step 7.3: Drop sequences too short to align full length

```bash
python3 07_solo1_homologues/scripts/7.3_filter_by_min_length.py \
  --input-fasta 07_solo1_homologues/data/results/7.2_solo1_clade_nonredundant.fasta \
  --output-fasta 07_solo1_homologues/data/results/7.3_solo1_clade_min800.fasta \
  --report-tsv 07_solo1_homologues/data/results/7.3_length_filter_report.tsv \
  --min-length 800
```

Outputs: [`7.3_solo1_clade_min800.fasta`](data/results/7.3_solo1_clade_min800.fasta),
86 sequences; 10 were shorter than 800 residues.

## Step 7.4: Align the full-length sequences

```bash
mafft --maxiterate 1000 --localpair \
  07_solo1_homologues/data/results/7.3_solo1_clade_min800.fasta \
  > 07_solo1_homologues/data/results/7.4_solo1_clade.mafft.fasta
```

## Step 7.5: Trim the alignment

```bash
clipkit 07_solo1_homologues/data/results/7.4_solo1_clade.mafft.fasta \
  -m gappy -g 0.9 \
  -o 07_solo1_homologues/data/results/7.5_solo1_clade.mafft.clipkit.fasta
```

## Step 7.6: Infer the tree

```bash
iqtree -s 07_solo1_homologues/data/results/7.5_solo1_clade.mafft.clipkit.fasta \
  -pre 07_solo1_homologues/data/results/7.6_solo1_clade \
  -st AA -m MFP -bb 1000 -nt 8
```

The model selected by BIC is `Q.BIRD+F+R6`.

Outputs: [`7.6_solo1_clade.treefile`](data/results/7.6_solo1_clade.treefile) with
ultrafast bootstrap support, the consensus tree
[`7.6_solo1_clade.contree`](data/results/7.6_solo1_clade.contree) and the run
report [`7.6_solo1_clade.iqtree`](data/results/7.6_solo1_clade.iqtree).

## Step 7.7: Measure the domain lengths

### Purpose

Measure the CC and NB-ARC regions of each homologue from its NLRtracker
annotation. The LRR region and the stretch from the end of the LRR to the C
terminus are measured on the alignment, between the columns that carry
positions 510 and 833 of the reference receptor `XP_020170543.1`, so that the
two are delimited the same way in every sequence.

### Script and command

```bash
python3 07_solo1_homologues/scripts/7.7_build_domain_length_table.py \
  --target-fasta 07_solo1_homologues/data/results/7.3_solo1_clade_min800.fasta \
  --alignment-fasta 07_solo1_homologues/data/results/7.5_solo1_clade.mafft.clipkit.fasta \
  --reference-seq-id Poales_Aegilops_tauschii_subsp._strangulata_XP_020170543.1_CNLO \
  --reference-lrr-start 510 --reference-lrr-end 833 \
  --seqid-assembly-tsv 07_solo1_homologues/data/inputs/7.0_seqid_to_assembly.tsv \
  --exon-boundaries-tsv 07_solo1_homologues/data/inputs/7.0_exon_boundaries.tsv \
  --nlrtracker-root 07_solo1_homologues/data/inputs/7.0_nlrtracker_domains \
  --output-table-tsv 07_solo1_homologues/data/results/7.7_domain_length_summary.tsv \
  --output-missing-tsv 07_solo1_homologues/data/results/7.7_domain_length_missing.tsv \
  --output-lrr-tail-fasta 07_solo1_homologues/work/7.7_lrr_tail.fasta \
  --output-itol-prefix 07_solo1_homologues/data/results/7.7_itol
```

### Outputs

- [`7.7_itol.domain_lengths_stacked.txt`](data/results/7.7_itol.domain_lengths_stacked.txt):
  the stacked CC, NB-ARC, LRR and C-terminal lengths shown in Fig. S5A.
- [`7.7_domain_length_summary.tsv`](data/results/7.7_domain_length_summary.tsv):
  the same measurements per sequence.

## Step 7.8: Assign each homologue to a grass clade

### Purpose

Resolve the species of each sequence against the NCBI taxonomy and assign it to
PACMAD or to one of the BOP subgroups, so that the branches of the tree can be
coloured by the grass lineage the receptor comes from.

### Command

```bash
python3 07_solo1_homologues/scripts/7.8_make_itol_bop_pacmad.py \
  --input-fasta 07_solo1_homologues/data/results/7.3_solo1_clade_min800.fasta \
  --output-species-taxonomy-tsv 07_solo1_homologues/data/results/7.8_species_taxonomy_taxonkit.tsv \
  --output-seqid-clade-tsv 07_solo1_homologues/data/results/7.8_seqid_to_bop_pacmad.tsv \
  --output-missing-tsv 07_solo1_homologues/data/results/7.8_taxonkit_missing.tsv \
  --output-itol-colorstrip 07_solo1_homologues/data/results/7.8_itol_bop_pacmad.colorstrip.txt
```

The script queries the NCBI taxonomy through `taxonkit`, which needs the
taxonomy files installed locally. Its outputs are distributed here so that the
rest of the workflow runs without them.

### Outputs

- [`7.8_itol_bop_pacmad.colorstrip.txt`](data/results/7.8_itol_bop_pacmad.colorstrip.txt):
  the branch colours of Fig. S5A.
- [`7.8_seqid_to_bop_pacmad.tsv`](data/results/7.8_seqid_to_bop_pacmad.tsv) and
  [`7.8_species_taxonomy_taxonkit.tsv`](data/results/7.8_species_taxonomy_taxonkit.tsv):
  the assignment of each sequence and the lineage it rests on.

## Step 7.9: Build the tip labels

```bash
python3 07_solo1_homologues/scripts/7.9_make_itol_labels.py \
  --input-fasta 07_solo1_homologues/data/results/7.3_solo1_clade_min800.fasta \
  --output-label-tsv 07_solo1_homologues/data/results/7.9_itol_labels.tsv \
  --output-itol-labels 07_solo1_homologues/data/results/7.9_itol_labels.txt
```

Outputs: [`7.9_itol_labels.txt`](data/results/7.9_itol_labels.txt), which
replaces each identifier with the species name and protein accession.

## Step 7.10: Build the alignment track

```bash
python3 07_solo1_homologues/scripts/7.10_make_itol_alignment.py \
  --input-alignment-fasta 07_solo1_homologues/data/results/7.5_solo1_clade.mafft.clipkit.fasta \
  --output-itol-alignment 07_solo1_homologues/data/results/7.10_itol_alignment.txt \
  --output-summary-tsv 07_solo1_homologues/data/results/7.10_itol_alignment_summary.tsv
```

Outputs: [`7.10_itol_alignment.txt`](data/results/7.10_itol_alignment.txt), the
alignment shown in Fig. S6.

## Step 7.11: Display the tree

Upload [`7.6_solo1_clade.treefile`](data/results/7.6_solo1_clade.treefile) to
[iTOL](https://itol.embl.de), then the labels of Step 7.9, the branch colours of
Step 7.8 and the domain-length track of Step 7.7, which together give Fig. S5A.
Uploading the alignment track of Step 7.10 instead gives Fig. S6. Fig. 2B shows
the six homologues that were tested, taken from the same tree.
