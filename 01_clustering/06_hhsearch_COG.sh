#add directory or directories containing 
ASGARDALID="data/cogs_v2 data/de_novo"
#Download COG database from ftp://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/hhsuite_dbs/COG_KOG.tar.gz
COGHHDB="data/COG_KOG/COG_KOG"
THREADS=15

find -name "*.aln" ${ASGARDALID} | awk 'NR>1' | parallel -j $THREADS "hhsearch -i {} -d $COGHHDB -glob -M 50 -blasttab stdout" > results/Asgard_db_COG_hhsearch.tsv
grep "^ \|Command" results/Asgard_db_COG_hhsearch.tsv | grep -v "No Hit" | sed 's/.*hhsuite\/Asgard_db\/aln\///' | sed 's/ -d .*//' | sed -E 's/^ [0-9]+ /\t/' | sed -e 's/^  [0-9] /\t/' | sed -E 's/ ...................... //' | sed -E 's/(100\.0) / \1 /' | sed -E 's/  +/ /'g  | sed
 's/ /\t/g' | awk -F'\t' '{if ($1 != "") target=$1; else print target,$0}' | sed  's/[0-9](/\t(/' > results/Asgard_db_COG_hhsearch_parse.tsv
grep -a "NAME  " ${COGHHDB}_hhm.ffdata | sed 's/NAME  //' | sed 's/ /\t/' > results/COG_KOG_function.tsv
