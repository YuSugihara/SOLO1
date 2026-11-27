# AlphaFold 3 modeling

This workflow converts the 693 NLR-IDs of the discovery set into the AlphaFold 3
predictions used for structural ranking. Each receptor is modeled as a
homohexamer in triplicate from its N terminus to the end of its NB-ARC domain.

## Inputs

- [`1.2_nlr_id_filtered.fasta`](../01_nlr_id_dataset/data/results/1.2_nlr_id_filtered.fasta)
  and
  [`1.2_nlr_id_filtered_domains.tsv`](../01_nlr_id_dataset/data/results/1.2_nlr_id_filtered_domains.tsv):
  the 693 NLR-IDs and their NLRtracker domain boundaries, from Step 1.2.

## Step 2.1: Extract the N-terminus-to-NB-ARC region for modeling

### Purpose

Extract the region from the N terminus to the end of the NB-ARC domain for
AlphaFold 3 modeling. Sequences containing an ambiguous residue (X) within that
region are excluded because the [AlphaFold Server](https://alphafoldserver.com)
rejects them. Each identifier
is formatted as `<accession>_6x<protein>_1-<end>`, with dots replaced by
underscores, where `6x` indicates hexameric submission.

### Script and command

```bash
python3 02_af3_modeling/scripts/2.1_extract_modelled_region.py \
  --input-fasta 01_nlr_id_dataset/data/results/1.2_nlr_id_filtered.fasta \
  --input-domains 01_nlr_id_dataset/data/results/1.2_nlr_id_filtered_domains.tsv \
  --output-fasta 02_af3_modeling/data/results/2.1_modelled_regions.fasta \
  --report-tsv 02_af3_modeling/data/results/2.1_modelled_regions_report.tsv
```

### Outputs

- [`2.1_modelled_regions.fasta`](data/results/2.1_modelled_regions.fasta): 689
  regions.
- [`2.1_modelled_regions_report.tsv`](data/results/2.1_modelled_regions_report.tsv):
  all 693 proteins with the modeled region and filter decision; 4 are excluded
  for containing an ambiguous residue.

## Step 2.2: Remove redundancy

### Purpose

Cluster the modeled regions at 95% sequence identity and keep one
representative per cluster so that closely related receptors are not modeled
repeatedly.

### Command

```bash
cd-hit \
  -i 02_af3_modeling/data/results/2.1_modelled_regions.fasta \
  -o 02_af3_modeling/data/results/2.2_nonredundant_regions.fasta \
  -c 0.95 -n 5
```

### Outputs

- [`2.2_nonredundant_regions.fasta`](data/results/2.2_nonredundant_regions.fasta):
  301 representatives.
- [`2.2_nonredundant_regions.fasta.clstr`](data/results/2.2_nonredundant_regions.fasta.clstr):
  cluster membership of all 689 regions.

## Step 2.3: Drop the representatives the AlphaFold Server cannot accept

### Purpose

Each representative is submitted as a hexamer together with 50 oleic acid
molecules (as a membrane proxy), so job size grows with region length. Five
representatives exceed the AlphaFold Server input size limit and cannot be
submitted. They are the five longest sequences, at 733–765 amino acids; the
longest accepted representative is 674 amino acids, so any threshold between
the two reproduces the selection. The script uses 700.

### Script and command

```bash
python3 02_af3_modeling/scripts/2.3_select_af3_targets.py \
  --input-fasta 02_af3_modeling/data/results/2.2_nonredundant_regions.fasta \
  --output-fasta 02_af3_modeling/data/results/2.3_af3_targets.fasta \
  --report-tsv 02_af3_modeling/data/results/2.3_af3_target_selection.tsv \
  --max-length 700
```

### Outputs

- [`2.3_af3_targets.fasta`](data/results/2.3_af3_targets.fasta): the 296 regions
  submitted for modeling.
- [`2.3_af3_target_selection.tsv`](data/results/2.3_af3_target_selection.tsv):
  all 301 representatives with their length and filter decision.

## Step 2.4: Write the AlphaFold Server job files

### Purpose

Write the 296 regions as AlphaFold Server batch job files. Each job contains
six copies of the protein chain, enables structural templates, and includes 50
oleic acid molecules (`CCD_OLA`) as a membrane proxy. Model seeds are left empty
so that the server assigns an independent seed to each submission; submitting the
same files once per replicate therefore yields independent predictions. The
server accepts at most 20 jobs per file.

### Script and command

```bash
python3 02_af3_modeling/scripts/2.4_generate_af3_jobs.py \
  --input-fasta 02_af3_modeling/data/results/2.3_af3_targets.fasta \
  --output-dir 02_af3_modeling/data/results/2.4_af3_jobs
```

### Outputs

- [`2.4_af3_jobs/`](data/results/2.4_af3_jobs): 296 jobs across 15 files.

## Step 2.5: Submit the jobs

The job files were submitted to the
[AlphaFold Server](https://alphafoldserver.com) three times, yielding three
independent predictions per receptor: 888 candidate predictions plus 3 for the
NbNRC2a benchmark submitted separately under identical conditions (891
predictions in total). The job requests returned by the server, which record the
seed assigned to each prediction, are kept as
[`2.5_submitted_job_requests.tar.gz`](data/results/2.5_submitted_job_requests.tar.gz).

The predicted structures themselves are too large to distribute through GitHub
and are deposited at Zenodo,
[doi:10.5281/zenodo.24540207](https://doi.org/10.5281/zenodo.24540207)
(`hexamer_models_rep1.tar.gz`–`hexamer_models_rep3.tar.gz`). The six structural
metrics computed from these predictions are provided in
[`03_structural_metrics`](../03_structural_metrics).
