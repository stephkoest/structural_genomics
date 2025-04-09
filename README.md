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
conda install -n foldseek -c biocondafoldseek::6.29e2557::pl5321hb365157_2
#colabfold has to be installed via pip
conda create -n colabfold
conda activate colabfold
pip install alphafold-colabfold::2.1.16
conda deactivate
```

### Other Software

Ensure you have the following software versions:

- MMseqs2 for ColabFold search: fd1837b600c57278bcfb2ac1ac7f024e458c0606

## Usage

Change into the respective subdirectory and run the scripts in ascending order.

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


### 04_prediction

1. [01_prepareQueries.sh](https://github.com/stephkoest/structural_genomics/blob/main/04_prediction/01_prepareQueries.sh)
   - Prepares queries for prediction by extracting and processing input files.
   - Activates the required Conda environment and runs `colabfold_split_msas` to generate MSA files.

2. [02_submit_prediction.sh](https://github.com/stephkoest/structural_genomics/blob/main/04_prediction/02_submit_prediction.sh)
   - Submits prediction jobs to a SLURM scheduler.
   - Activates the required Conda environment and runs `colabfold_batch` for predictions.

3. [03_sum_predictions.sh](https://github.com/stephkoest/structural_genomics/blob/main/04_prediction/03_sum_predictions.sh)
   - Summarizes prediction results by extracting and processing prediction logs and PDB files.
   - Calculates average pLDDT scores for predicted structures and selects the best structures.

4. [04_select_leftovers.sh](https://github.com/stephkoest/structural_genomics/blob/main/04_prediction/04_select_leftovers.sh)
   - Identifies and selects sequences that were not successfully predicted.
   - Prepares a list of missing sequences for further processing or prediction attempts.


### 05_structure_search

1. [00_prepare_db.sh](https://github.com/stephkoest/structural_genomics/blob/main/05_structure_search/00_prepare_db.sh)
   - Creates necessary directories (data, results, log).
   - Activates the foldseek Conda environment.
   - Prepares a database from PDB files located in ../04_prediction/results/Asgard_best_structures/.
     
2. [01_submit_foldseek.sh](https://github.com/stephkoest/structural_genomics/blob/main/05_structure_search/01_submit_foldseek.sh)
   - Configures a SLURM job for structure search.
   - Activates the foldseek Conda environment.
   - Iterates over multiple databases to perform foldseek searches and converts alignments to m8 format.
   - Outputs results to the results directory.

3. [03_cluster_CIPS.sh](https://github.com/stephkoest/structural_genomics/blob/main/05_structure_search/03_cluster_CIPS.sh)
   - Activates the foldseek environment (commented out in the script).
   - Performs foldseek search on the Asgard database against itself.
   - Generates TSV files for alignments and clusters, storing results in the results directory.

### 06_gene_tree_annotation

1. [generate_input_data_iTol.py](https://github.com/stephkoest/structural_genomics/blob/main/06_gene_tree_annotation/generate_input_data_iTol.py)
   - Generates datasets to annotate and visualize single gene trees in iTOL
   - See its [README.md](https://github.com/stephkoest/structural_genomics/blob/main/06_gene_tree_annotation/README.md) for further details

## Contributors

- [Stephan Köstlbacher](https://github.com/stephkoest)
- [Jolien van Hooff](https://github.com/jolienvanhooff)

## License

This project is licensed under the MIT License.
