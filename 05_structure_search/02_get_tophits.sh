#!/bin/bash
#get non-overlapping best hits per query structure
for file in results/*m8
do
	OUTFILE=results/$(basename $file .m8)_top.m8
	Rscript bin/best_hits.R -i $file -o $OUTFILE
done
