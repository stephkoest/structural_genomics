#!/bin/bash
conda activate mmseqs2_v14
mmseqs createdb data/asCOG_seq_domains.faa asCOG_seq_domains_db
mkdir -p log
sbatch bin/runMmseqsProfile.sh data/Asgard_DB asCOG_seq_domains_db
