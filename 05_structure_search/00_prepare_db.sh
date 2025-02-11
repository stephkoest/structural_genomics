#!/bin/bash
mkdir -p data results
#load foldseek
conda activate foldseek

foldseek createdb ../04_prediction/results/Asgard_best_structures/*.pdb data/Asgard_db

