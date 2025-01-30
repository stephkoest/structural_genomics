#!/bin/bash
conda activate mmseqs2_v14
mmseqs createdb data/Asgard_DB.faa data/Asgard_DB
mmseqs createdb data/asCOG_seq_domains.faa data/asCOG_seq_domains_db
grep ">" data/Asgard_DB.faa | sed 's/>//' | cut -d" " -f1 > data/Asgard_DB.txt
mkdir -p log
sbatch bin/runMmseqsProfile.sh data/Asgard_DB asCOG_seq_domains_db
