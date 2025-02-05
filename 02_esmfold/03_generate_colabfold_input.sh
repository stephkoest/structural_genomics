#!/bin/bash
OUTD=../03_genetic_search/data/Asgard_DB
mkdir -p $(dirname $OUTD) $OUTD
grep -A1 -Ff <(awk '$2<80 {print $1}' results/Asgard_DB_de_novo_representatives_plDDT.tsv | sed 's/\.pdb//' | rev | cut -f1 -d"/" | rev) data/Asgard_DB_de_novo_representatives.faa | grep -v "^--$" > $OUTD/Asgard_DB_remaining_representatives.faa
grep -A1 -Ff <(awk '$2<80 {print $1}' results/Asgard_DB_asCOG_representatives_plDDT.tsv | sed 's/\.pdb//' | rev | cut -f1 -d"/" | rev) data/Asgard_DB_asCOG_representatives.faa | grep -v "^--$" >> $OUTD/Asgard_DB_remaining_representatives.faa
