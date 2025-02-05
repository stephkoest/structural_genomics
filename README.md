# Structural Genomics

Collection of scripts used for analyses in "Structure-based inference of eukaryotic complexity in Asgard archaea."

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Setup](#setup)
- [Usage](#usage)
- [Subdirectory Summaries](#subdirectory-summaries)
- [Contributors](#contributors)
- [License](#license)

## Introduction

This repository contains a collection of scripts and tools used for various analyses as described in the paper "Structure-based inference of eukaryotic complexity in Asgard archaea."

## Features

- MMseqs2 integration for genetic search.
- Adapted ColabFold genetic search script with support for external databases.
- Various utilities for data analysis and visualization.

## Setup

### Conda Environments

To set up the necessary conda environments, run:

```bash
conda install -n mmseqs2_v14 -c bioconda mmseqs2::14.7e284
```

### Other Software

Ensure you have the following software versions:

- MMseqs2 for ColabFold search: fd1837b600c57278bcfb2ac1ac7f024e458c0606

## Usage

to be done

## Subdirectory Summaries

### 01_clustering

1. [01_submitMmseqs.sh](https://github.com/stephkoest/structural_genomics/blob/main/01_clustering/01_submitMmseqs.sh)
   - This script submits MMseqs2 jobs for clustering sequences.

2. [02_cluster_map.R](https://github.com/stephkoest/structural_genomics/blob/main/01_clustering/02_cluster_map.R)
   - An R script for mapping the clusters obtained from MMseqs2.

3. [03_extractAlignDomains.sh](https://github.com/stephkoest/structural_genomics/blob/main/01_clustering/03_extractAlignDomains.sh)
   - This script extracts and aligns domain sequences from the clustered data.

4. [04_cluster_leftover.sh](https://github.com/stephkoest/structural_genomics/blob/main/01_clustering/04_cluster_leftover.sh)
   - Clusters any leftover sequences that were not included in the initial clustering.

5. [05_select_cluster_reps.R](https://github.com/stephkoest/structural_genomics/blob/main/01_clustering/05_select_cluster_reps.R)
   - An R script to select representative sequences from each cluster.

6. [06_hhsearch_COG.sh](https://github.com/stephkoest/structural_genomics/blob/main/01_clustering/06_hhsearch_COG.sh)
   - This script performs an HHsearch against the COG database to annotate the clusters.

### 02_esmfold

1. [01_submit_ESMfold.sh](https://github.com/stephkoest/structural_genomics/blob/main/02_esmfold/01_submit_ESMfold.sh)
   - This script submits ESMfold jobs to a cluster. It processes input protein sequence files and splits them into smaller chunks, then submits these chunks as separate jobs for ESMfold prediction.

2. [02_sum_ESMfold.sh](https://github.com/stephkoest/structural_genomics/blob/main/02_esmfold/02_sum_ESMfold.sh)
   - This script summarizes the results of ESMfold predictions by calculating the average plDDT score for each protein structure.

3. [03_generate_colabfold_input.sh](https://github.com/stephkoest/structural_genomics/blob/main/02_esmfold/03_generate_colabfold_input.sh)
   - This script generates input files for ColabFold based on the summarized ESMfold results. It filters sequences with an average plDDT score below a certain threshold and prepares them for further analysis.

### 03_genetic_search

1. [01_create_chunks.sh](https://github.com/stephkoest/structural_genomics/blob/main/03_genetic_search/01_create_chunks.sh)
   - This script processes input protein sequence files, categorizes them based on sequence length, and creates compressed chunks for further analysis.

2. [02_submitSearch.sh](https://github.com/stephkoest/structural_genomics/blob/main/03_genetic_search/02_submitSearch.sh)
   - This script submits genetic search jobs to a cluster, utilizing MMseqs2 for sequence alignment and searching against predefined databases. It also manages result packaging and compression.

## Contributors

- [Stephan Koestlbacher](https://github.com/stephkoest)

## License

This project is licensed under the MIT License.
