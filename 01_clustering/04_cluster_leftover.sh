#!/bin/bash
LEFT=results/Asgard_DB_asCOG_left_partial_proteins.tsv
UNMAPPED=results/Asgard_DB_asCOG_unmapped.txt
INFASTA=data/Asgard_DB_00v03_01v04_shortHeader.faa
OUTFASTA=results/Asgard_DB_asCOG_left.faa
CLUSTBASE=results/Asgard_DB_left_20_blast_c50

### Extract reps
ASCOGREPS=results/Asgard_DB_asCOG_representatives.tsv
REPSEQS=results/Asgard_DB_asCOG_representatives.faa
grep -A1 -wFf <(cut -f2 $ASCOGREPS | sed -E 's/__[0-9]+-[0-9]+//' | grep -v "^NA$" | sort -u) $INFASTA \
	| grep -v "^--$" > $REPSEQS
###

grep -A1 -Ff <(awk 'NR>1{print $1}' $LEFT) $INFASTA \
    | grep -v "^--$" \
    | python bin/extractMultiDomain.py <(awk 'NR>1{print $1","$2"-"$3}' $LEFT) > ${OUTFASTA}

grep -A1 -Ff <(awk '{print $1}' $UNMAPPED) $INFASTA \
    | grep -v "^--$" >> ${OUTFASTA}

#cluster leftovers
mmseqs easy-cluster ${OUTFASTA} $CLUSTBASE tmp \
	--cluster-reassign \
	--min-seq-id 0.2 \
	--split-memory-limit 110G  \
	--cov-mode 0 \
	--cluster-mode 1 \
	-e 0.001 -c 0.5 --single-step-clustering 1

mkdir -p data/de_novo
#Extract cluster representatives
parallel -j 38 "grep -A1 -Ff <(grep {} ${CLUSTBASE}_cluster.tsv | cut -f2)  $OUTFASTA | grep -v \"^--$\" > data/de_novo/{}.faa" ::: $(awk '{print $1}' ${CLUSTBASE}_cluster.tsv | sort | uniq -c | awk '$1>4 {print $2}')

conda activate famsa2
parallel -j 20 "echo working on {} && famsa -t 1 -refine_mode on {} data/de_novo/{/.}.aln" ::: data/de_novo/*.faa

mkdir -p data/de_novo_sto
conda activate seqmagick

parallel -j 20 "seqmagick convert --input-format fasta --output-format stockholm {} data/de_novo_sto/{/.}.sto" ::: data/de_novo/*.aln

conda activate mmseqs2
for file in data/de_novo_sto/*.sto
do
        echo "# STOCKHOLM 1.0"
        echo "#=GF AC $(basename $file .sto)"
        awk 'NR>1' $file | grep -v "#=GS"
        echo "//"
done | gzip > data/de_novo_alpaca_msa.gz

rm -r stockholm

mmseqs convertmsa data/de_novo_alpaca_msa.gz data/de_novo_alpaca_msa_db
mmseqs msa2profile data/de_novo_alpaca_msa_db data/de_novo_alpaca_profile_db --match-mode 1
cat data/de_novo/*.faa > data/de_novo_sto_alpaca.faa
mmseqs createdb data/de_novo_sto_alpaca.faa data/de_novo_alpaca_db
mmseqs search data/de_novo_alpaca_db data/de_novo_alpaca_profile_db data/de_novo_alpaca_res tmp -e 0.001 -s 9 --remove-tmp-files 1 --threads 38 --split-memory-limit 80G
mmseqs convertalis data/de_novo_alpaca_db data/de_novo_alpaca_profile_db data/de_novo_alpaca_res data/de_novo_alpaca_res.m8 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
