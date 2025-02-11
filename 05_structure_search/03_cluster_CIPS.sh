#!/bin/bash
#source /local/one/software/miniconda3/bin/activate
#conda activate foldseek
DB=data/Asgard_db

foldseek search $DB $DB aln tmpFolder -a
foldseek createtsv $DB $DB aln results/Asgard_db_aln.tsv
foldseek clust $DB aln clu
foldseek createtsv $DB $DB clu results/Asgard_db_cluster.tsv
