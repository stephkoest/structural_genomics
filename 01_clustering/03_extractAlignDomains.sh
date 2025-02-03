#!/bin/bash
source /local/one/software/miniconda3/bin/activate
COGASS=results/Asgard_DB_asCOG_safe_v2.tsv 
INFASTA=data/Asgard_DB_00v03_01v04.faa
OUTFASTA=data/$(basename $INFASTA .faa)_shortHeader.faa
sed 's/ .*//' $INFASTA \
    | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' \
    | sed 's/\t/\n/' \
    > $OUTFASTA
INFASTA=$OUTFASTA
mkdir -p data/cogs_v2

parallel -j 38 "bash bin/grepDomains.sh {} $COGASS $OUTFASTA data/cogs_v2" ::: $(awk 'NR>1 {print $1}' $COGASS | sort -u )
conda activate famsa2
parallel -j 38 "echo working on {} && famsa -t 1 -refine_mode on {} data/cogs_v2/{/.}.aln" ::: data/cogs_v2/*.faa

mkdir -p data/cogs_v2_sto
conda activate seqmagick

parallel -j 38 "seqmagick convert --input-format fasta --output-format stockholm {} data/cogs_v2_sto/{/.}.sto" ::: data/cogs_v2/*.aln

conda activate mmseqs2
for file in data/cogs_v2_sto/*.sto
do
        echo "# STOCKHOLM 1.0"
        echo "#=GF AC $(basename $file .sto)"
        awk 'NR>1' $file | grep -v "#=GS"
        echo "//"
done | gzip > data/asCOG_alpaca_v2_msa.gz

rm -r data/cogs_v2_sto

mmseqs convertmsa data/asCOG_alpaca_v2_msa.gz data/asCOG_alpaca_v2_msa_db
mmseqs msa2profile data/asCOG_alpaca_v2_msa_db data/asCOG_alpaca_v2_profile_db --match-mode 1
cat data/cogs_v2/cog.0*.faa > data/cogs_v2_alpaca.faa
mmseqs createdb data/cogs_v2_alpaca.faa data/cogs_v2_alpaca_db
mmseqs search data/cogs_v2_alpaca_db data/asCOG_alpaca_v2_profile_db data/cogs_v2_alpaca_res tmp -e 0.001 -s 9 --remove-tmp-files 1 --threads 38 --split-memory-limit 80G
mmseqs convertalis data/cogs_v2_alpaca_db data/asCOG_alpaca_v2_profile_db data/cogs_v2_alpaca_res data/cogs_v2_alpaca_res.m8 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
