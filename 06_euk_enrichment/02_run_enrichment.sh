#!/bin/bash
mkdir -p results
Rscript run_euk_enrichment.R \
  --aln ../05_structure_search/results/Asgard_UniProt50_Alphafold_aln.m8 \
  --map data/uniref50.sub.long.tsv \
  --tax data/uniref50.sub.tsv \
  --db data/accessionTaxa.sql \
  --out results/euk_enrichment_summary
