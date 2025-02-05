#!/bin/bash
find results/chunks/ -name "*asCOG*.pdb" -exec awk -v prot=$(basename {} .pdb) 'NR>1{sum += $11}END{print prot"\t"sum/NR}' {} \; > results/Asgard_DB_asCOG_representatives_plDDT.tsv
find results/chunks/ -name "*de_novo*.pdb" -exec awk -v prot=$(basename {} .pdb) 'NR>1{sum += $11}END{print prot"\t"sum/NR}' {} \; > results/Asgard_DB_de_novo_representatives_plDDT.tsv
