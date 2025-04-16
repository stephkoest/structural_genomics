INFASTA=data/Asgard_DB_00v03_01v04_shortHeader.faa
REPSEQS=results/Asgard_DB_de_novo_representatives.faa
#for the de novo clusters
grep -A1 -wFf <(cut -f1 results/Asgard_DB_de_novo_representatives.tsv | sed -E 's/__[0-9]+-[0-9]+//' | grep -v "^NA$" | sort -u) $INFASTA \
        | grep -v "^--$" > $REPSEQS
#for the ascogs
REPSEQS=results/Asgard_DB_asCOG_representatives.faa
grep -A1 -wFf <(cut -f2 results/Asgard_DB_asCOG_representatives.tsv | sed -E 's/__[0-9]+-[0-9]+//' | grep -v "^NA$" | sort -u) $INFASTA \
        | grep -v "^--$" > $REPSEQ
