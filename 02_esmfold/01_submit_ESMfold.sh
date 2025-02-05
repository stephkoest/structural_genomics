#!/bin/bash
IND=../01_clustering/data
DATAD=data/chunks
RESD=results/chunks
mkdir -p $DATAD results $RESD
split -d -l 2000 ${IND}/Asgard_DB_de_novo_representatives.faa ${DATAD}/Asgard_DB_de_novo_representatives
split -d -l 2000 ${IND}/Asgard_DB_asCOG_representatives.faa ${DATAD}/Asgard_DB_asCOG_representatives

for file in ${DATAD}/Asgard_DB_*
do
	sbatch bin/ESMfold.sh $file ${RESD}/$(basename $file)
done
