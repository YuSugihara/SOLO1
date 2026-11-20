# NLR-ID discovery set

This workflow collects the coiled-coil NLRs that carry an integrated domain (ID)
from the NLRs annotated with NLRtracker in 180 RefSeq reference plant proteomes,
and removes the candidates whose extra domain is an annotation artefact. The
proteomes are the RefSeq annotations held at NCBI; the NLRtracker annotations
were produced and deposited by
[Toghani and Kamoun (2024)](https://doi.org/10.5281/zenodo.13627395).

The 180 assemblies, their species and order, proteome size and NLR count are
listed in [`refseq_assemblies.tsv`](data/inputs/refseq_assemblies.tsv). Each
assembly corresponds to a different species, so the 180 proteomes represent 180
species. Across them NLRtracker annotates 113,684 NLRs.

## Inputs

The per-assembly NLRtracker output directories are not redistributed here.
Download the deposited dataset from
[doi:10.5281/zenodo.13627395](https://doi.org/10.5281/zenodo.13627395) and point
`NLRTRACKER_DIR` at the directory that holds them, each named
`nlrtracker_<assembly accession>_<suffix>`.

```bash
NLRTRACKER_DIR=/path/to/RefSeq_NLRtracker
```

- [`refseq_assemblies.tsv`](data/inputs/refseq_assemblies.tsv): the 180
  assembly accessions with species, order, proteome size and NLR count, taken
  from Data S2 of the same deposit.
- [`refseq_nlrtracker_tables.tsv.gz`](data/inputs/refseq_nlrtracker_tables.tsv.gz):
  the NLRtracker classification tables of all 180 assemblies concatenated, with
  the assembly accession added as the first column; 128,425 rows. Step 1.1 reads
  the original per-assembly directories rather than this file, which is provided
  so that the protein counts can be checked without downloading the archive:

```bash
gzip -cd 01_nlr_id_dataset/data/inputs/refseq_nlrtracker_tables.tsv.gz \
  | awk -F'\t' 'NR>1 && $3=="NLR" {n++} NR>1 && $3=="NLR" && $6=="CNLO" {c++} END {print n, c}'
# 113684 938
```

## Step 1.1: Extract CC-NLRs that carry an extra domain

### Purpose

Retain every protein that NLRtracker classifies with `Status = NLR` and
simplified domain architecture `CNLO`, that is, a coiled-coil NLR carrying at
least one domain in addition to CC, NB-ARC and LRR. The assembly accession is
prefixed to every sequence identifier and to every annotation row, so that each
record stays traceable to its genome. Assemblies and proteins are processed in
sorted order, so the output does not depend on the order in which the
filesystem lists the input directories.

### Script and command

```bash
python3 01_nlr_id_dataset/scripts/1.1_extract_nlr_id_candidates.py \
  --nlrtracker-dir "$NLRTRACKER_DIR" \
  --output-dir 01_nlr_id_dataset/data/results \
  --report-tsv 01_nlr_id_dataset/data/results/1.1_nlr_id_candidates_per_assembly.tsv
```

### Outputs

- [`1.1_nlr_id_candidates.fasta`](data/results/1.1_nlr_id_candidates.fasta): 938
  full-length protein sequences.
- [`1.1_nlr_id_candidates_nbarc.fasta`](data/results/1.1_nlr_id_candidates_nbarc.fasta):
  938 NB-ARC sequences.
- [`1.1_nlr_id_candidates.gff`](data/results/1.1_nlr_id_candidates.gff):
  InterProScan signatures of the 938 candidates. The GFF seqid is
  `<assembly accession>|<protein identifier>`.
- [`1.1_nlr_id_candidates_domains.tsv`](data/results/1.1_nlr_id_candidates_domains.tsv):
  NLRtracker domain table of the 938 candidates, with the assembly accession
  added as a first column. `type` is `REGION` or `DOMAIN` and `description`
  names the region or domain, so the NB-ARC region of a protein is the row with
  `type = REGION` and `description = NBARC`.
- [`1.1_nlr_id_candidates_per_assembly.tsv`](data/results/1.1_nlr_id_candidates_per_assembly.tsv):
  NLR and candidate counts per assembly; 180 rows. 87 assemblies contribute at
  least one candidate.

## Step 1.2: Remove annotation artefacts

### Purpose

Two classes of candidate are removed:

- **Transposon.** The protein carries a domain annotated as `TRANSPOSON`.
- **LRR-internal.** Every extra domain of the protein falls between the start of
  the first and the end of the last annotated LRR.

A protein is retained when it carries at least one extra domain outside the LRR
region and no `TRANSPOSON` domain. The report records every candidate once, with
its decision, the LRR span and the coordinates of each extra domain.

### Script and command

```bash
python3 01_nlr_id_dataset/scripts/1.2_filter_nlr_id_candidates.py \
  --input-dir 01_nlr_id_dataset/data/results \
  --output-dir 01_nlr_id_dataset/data/results \
  --report-tsv 01_nlr_id_dataset/data/results/1.2_nlr_id_filter_report.tsv
```

### Outputs

- [`1.2_nlr_id_filtered.fasta`](data/results/1.2_nlr_id_filtered.fasta): 693
  full-length NLR-ID sequences.
- [`1.2_nlr_id_filtered_nbarc.fasta`](data/results/1.2_nlr_id_filtered_nbarc.fasta):
  693 NB-ARC sequences.
- [`1.2_nlr_id_filtered.gff`](data/results/1.2_nlr_id_filtered.gff) and
  [`1.2_nlr_id_filtered_domains.tsv`](data/results/1.2_nlr_id_filtered_domains.tsv):
  signatures and domain boundaries of the 693 NLR-IDs. The domain table defines
  the modelled region used in the next analysis.
- [`1.2_nlr_id_filter_report.tsv`](data/results/1.2_nlr_id_filter_report.tsv): the
  decision for each of the 938 candidates; 938 rows.

### Counts

| Set | Proteins |
| --- | --- |
| NLRs annotated across the 180 proteomes | 113,684 |
| CC-NLRs carrying an extra domain (`CNLO`) | 938 |
| Removed: `TRANSPOSON` domain | 194 |
| Removed: every extra domain inside the LRR region | 51 |
| **Retained NLR-IDs** | **693** |
