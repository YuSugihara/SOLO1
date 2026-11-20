# AlphaFold 3 modelling

This workflow turns the 693 NLR-IDs of the discovery set into the AlphaFold 3
predictions that the structural ranking is built on. Each receptor is modelled
from its N terminus to the end of its NB-ARC domain, as a hexamer, in
triplicate.

## Inputs

- [`1.2_nlr_id_filtered.fasta`](../01_nlr_id_dataset/data/results/1.2_nlr_id_filtered.fasta)
  and
  [`1.2_nlr_id_filtered_domains.tsv`](../01_nlr_id_dataset/data/results/1.2_nlr_id_filtered_domains.tsv):
  the 693 NLR-IDs and their NLRtracker domain boundaries, from Step 1.2.

## Step 2.1: Cut each NLR-ID down to the modelled region

### Purpose

Take the region from the N terminus to the end of the NB-ARC domain, which is
the region submitted for modelling. Drop sequences that contain an ambiguous
residue (X) within that region, because the AlphaFold Server rejects them.
Rewrite each identifier as `<accession>_6x<protein>_1-<end>`, with dots
replaced by underscores, where `6x` records that the region is submitted as a
hexamer.

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
  all 693 proteins with the modelled region and the decision; 4 are dropped for
  containing an ambiguous residue.

## Step 2.2: Remove redundancy

### Purpose

Cluster the modelled regions at 95% sequence identity and keep one
representative per cluster, so that closely related receptors are not modelled
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
  the cluster membership of all 689 regions.

## Step 2.3: Drop the representatives the AlphaFold Server cannot accept

### Purpose

Each representative is submitted as a hexamer together with 50 copies of a
ligand, so the size of a job grows with the length of the region. Five
representatives exceed the size limit of the AlphaFold Server and cannot be
submitted. They are the five longest, at 733 to 765 residues; the longest
representative that is accepted is 674 residues, so any threshold between the
two reproduces the selection. The script uses 700.

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
  submitted for modelling.
- [`2.3_af3_target_selection.tsv`](data/results/2.3_af3_target_selection.tsv):
  all 301 representatives with their length and the decision.

## Step 2.4: Write the AlphaFold Server job files

### Purpose

Write the 296 regions as AlphaFold Server batch job files. Every job holds six
copies of the protein chain, uses structure templates, and includes 50 copies
of the oleic acid ligand `CCD_OLA`. Model seeds are left empty, so the server
assigns an independent seed to each submission; submitting the same files once
per replicate therefore gives independent predictions. The server accepts at
most 20 jobs per file.

### Script and command

```bash
python3 02_af3_modeling/scripts/2.4_generate_af3_jobs.py \
  --input-fasta 02_af3_modeling/data/results/2.3_af3_targets.fasta \
  --output-dir 02_af3_modeling/data/results/2.4_af3_jobs
```

### Outputs

- [`2.4_af3_jobs/`](data/results/2.4_af3_jobs): 296 jobs across 15 files.

## Step 2.5: Submit the jobs

The job files were submitted to the AlphaFold Server three times, giving three
independent predictions per receptor, 888 in total. The job requests returned
by the server, which record the seed assigned to each prediction, are kept as
[`2.5_submitted_job_requests.tar.gz`](data/results/2.5_submitted_job_requests.tar.gz).

The predicted structures themselves are too large to distribute through GitHub
and are deposited separately. The confidence summaries they are scored on are
held in the next analysis.
