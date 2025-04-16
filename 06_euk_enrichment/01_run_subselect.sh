#!/bin/bash
zgrep "entry id=\"UniRef50_\|\"common taxon\|property type=\"UniProtKB " data/uniref50.xml.gz > uniref50.sub.xml
sed 's/<entry id="//' uniref50.sub.xml \
	| sed 's/" updated="2022-12-14">//' \
	| sed 's/.*UniProtKB accession" value="/,/' \
	| sed 's/.*value="//' \
	| sed 's/".*//' \
	| sed 's/$/##/' \
	| sed 's/UniRef50_/\n/' \
	| awk 'BEGIN{RS="##\n"; ORS="\t"} {print }' \
	| sed -E 's/\t([0-9]+)\t,/\t\1\t/g' \
	| sed 's/\t,/,/g' \
	| grep -v "^$" \
	| awk -F"\t" 'OFS="\t" {if ($4 == "") print $0,$1; else print $0}' > data/uniref50.sub.tsv

rm uniref50.sub.xml
cat data/uniref50.sub.tsv | python bin/makeLong.py  > data/uniref50.sub.long.tsv
