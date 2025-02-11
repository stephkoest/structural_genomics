#!/bin/bash
awk '{print $1}' ../02_esmfold//Asgard_DB_*_representatives_plDDT.tsv results/Asgard_COLAB_pLDDT.tsv  |  rev | cut -d"/" -f1 | rev | sed 's/\..*//' | sed 's/_unrelaxed.*/
/' | sort -u > results/Asgard_DB_done.txt
grep -h ">" ../01_clustering/data/*_representatives.faa  | grep -vFf results/Asgard_DB_done.txt | sed 's/>//' > results/Asgard_DB_missing.txt
grep -A1 -whFf results/Asgard_DB_missing.txt ../01_clustering/data/*_representatives.faa | grep -v "^--$" > ../02_esmfold/data/Asgard_DB/Asgard_DB_remaining_representatives.faa
