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

7. [bin](https://github.com/stephkoest/structural_genomics/tree/main/01_clustering/bin)
   - A directory that potentially contains additional binaries or scripts used in the clustering process.

## Contributors

- [Stephan Koestlbacher](https://github.com/stephkoest)

## License

This project is licensed under the MIT License.
