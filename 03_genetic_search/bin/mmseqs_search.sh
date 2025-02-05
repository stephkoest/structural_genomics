#!/bin/bash
MMSEQS=$1
INFILE=$2
DB1=$3
DB2=$4
DB3=$5
CPU=$6

USE_ENV="1"
USE_EXT="1"

FILTER="1"
INBASE=$(basename $INFILE .faa)
BASE=${TMPDIR}/search
OUTD=$(dirname $INFILE)
INDB=${BASE}/${INBASE}_db
mkdir -p $BASE $OUTD

DB1ALI=${OUTD}/${INBASE}_uniref.a3m
DB2ALI=${OUTD}/${INBASE}_bfd.mgnify30.metaeuk30.smag30.a3m
DB3ALI=${OUTD}/${INBASE}_alpaca.a3m
MERGEALI=${OUTD}/${INBASE}_merged.a3m
MERGEPROF=${OUTD}/${INBASE}_merged_prof
FINALALI=${OUTD}/${INBASE}_merged_alpaca.a3m


SENSITIVITY=8
EXPAND_EVAL=0.1
DIFF=3000
ALIGN_EVAL=10
MAX_ACCEPT=100000
QSC=0.8

SEARCH_PARAM="--num-iterations 3 --db-load-mode 3 -a -s ${SENSITIVITY} -e 0.1 --max-seqs 10000 --threads ${CPU}"
FILTER_PARAM="--filter-msa ${FILTER} --filter-min-enable 1000 --diff ${DIFF} --qid 0.0,0.2,0.4,0.6,0.8,1.0 --qsc 0 --max-seq-id 0.95 --threads ${CPU}"
EXPAND_PARAM="--expansion-mode 0 -e ${EXPAND_EVAL} --expand-filter-clusters ${FILTER} --threads ${CPU} --max-seq-id 0.95"
OUTDB=${BASE}/${INBASE}_out
LTMPDIR=${TMPDIR}/${INBASE}


"$MMSEQS" createdb "$INFILE" "$INDB"
"$MMSEQS" search "$INDB" "$DB1" "$OUTDB" $LTMPDIR $SEARCH_PARAM
"$MMSEQS" expandaln "$INDB" "${DB1}.idx" "$OUTDB" "${DB1}.idx" "${OUTDB}_exp" ${EXPAND_PARAM}
"$MMSEQS" mvdb "${LTMPDIR}/latest/profile_1" "${BASE}/${INBASE}_prof_res"
"$MMSEQS" lndb "${INDB}_h" "${BASE}/${INBASE}_prof_res_h"
"$MMSEQS" align "${BASE}/${INBASE}_prof_res" "${DB1}.idx" "${OUTDB}_exp" "${OUTDB}_exp_realign" --db-load-mode 2 -e ${ALIGN_EVAL} --max-accept ${MAX_ACCEPT} --alt-ali 10 -a --threads ${CPU}
"$MMSEQS" filterresult "$INDB" "${DB1}.idx" "${OUTDB}_exp_realign" "${OUTDB}_exp_realign_filter"  --qid 0 --qsc $QSC --diff 0 --max-seq-id 1.0 --filter-min-enable 100 --db-load-mode 2 --threads ${CPU}
"$MMSEQS" result2msa "$INDB" "${DB1}.idx" "${OUTDB}_exp_realign_filter" $DB1ALI --msa-format-mode 6 ${FILTER_PARAM} --db-load-mode 2
"$MMSEQS" convertalis "${BASE}/${INBASE}_prof_res" $DB1 "${OUTDB}_exp_realign_filter" "${OUTDB}.m8" --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar 
"$MMSEQS" rmdb "${OUTDB}_exp_realign"
"$MMSEQS" rmdb "${OUTDB}_exp"
"$MMSEQS" rmdb "${OUTDB}"
"$MMSEQS" rmdb "${OUTDB}_exp_realign_filter"

if [ "${USE_ENV}" = "1" ]; then
  OUTDB=${BASE}/${INBASE}_out_env
  "$MMSEQS" search "${BASE}/${INBASE}_prof_res" "${DB2}" "${BASE}/${INBASE}_res_env" "${LTMPDIR}" $SEARCH_PARAM
  "$MMSEQS" expandaln "${BASE}/${INBASE}_prof_res" "${DB2}.idx" "${BASE}/${INBASE}_res_env" "${DB2}.idx" "${BASE}/${INBASE}_res_env_exp" -e ${EXPAND_EVAL} --expansion-mode 0 --db-load-mode 2
  "$MMSEQS" mvdb "${LTMPDIR}/latest/profile_1" "${BASE}/${INBASE}_prof_res_env"
  "$MMSEQS" lndb "${INDB}_h" "${BASE}/${INBASE}_prof_res_env_h"
  "$MMSEQS" align "${BASE}/${INBASE}_prof_res_env" "${DB2}.idx" "${BASE}/${INBASE}_res_env_exp" "${BASE}/${INBASE}_res_env_exp_realign" --db-load-mode 2 -e ${ALIGN_EVAL} --max-accept ${MAX_ACCEPT} --alt-ali 10 -a
  "$MMSEQS" filterresult "${INDB}" "${DB2}.idx" "${BASE}/${INBASE}_res_env_exp_realign" "${BASE}/${INBASE}_res_env_exp_realign_filter" --db-load-mode 2 --qid 0 --qsc $QSC --diff 0 --max-seq-id 1.0 --filter-min-enable 100
  "$MMSEQS" result2msa "${INDB}" "${DB2}.idx" "${BASE}/${INBASE}_res_env_exp_realign_filter" $DB2ALI --msa-format-mode 6 --db-load-mode 2 ${FILTER_PARAM}
  "$MMSEQS" convertalis "${BASE}/${INBASE}_prof_res_env" $DB2 "${BASE}/${INBASE}_res_env_exp_realign_filter" "${OUTDB}.m8" --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar --db-load-mode 2
  "$MMSEQS" mergedbs $DB2ALI $MERGEALI $DB1ALI $DB2ALI
  "$MMSEQS" rmdb "${BASE}/${INBASE}_res_env_exp_realign_filter"
  "$MMSEQS" rmdb "${BASE}/${INBASE}_res_env_exp_realign"
  "$MMSEQS" rmdb "${BASE}/${INBASE}_res_env_exp"
  "$MMSEQS" rmdb "${BASE}/${INBASE}_res_env"
fi

if [ "${USE_EXT}" = "1" ]; then
  OUTDB=${BASE}/${INBASE}_out_env_ext
  "$MMSEQS" search "${BASE}/${INBASE}_prof_res" "${DB3}" "${BASE}/${INBASE}_res_env_ext" "${LTMPDIR}" $SEARCH_PARAM
  "$MMSEQS" mvdb "${LTMPDIR}/latest/profile_1" "${BASE}/${INBASE}_prof_res_env_ext"
  "$MMSEQS" lndb "${INDB}_h" "${BASE}/${INBASE}_prof_res_env_ext_h"
  "$MMSEQS" filterresult "${INDB}" "${DB3}.idx" "${BASE}/${INBASE}_res_env_ext" "${BASE}/${INBASE}_res_env_ext_filter" --db-load-mode 2 --qid 0 --qsc $QSC --diff 0 --max-seq-id 1.0 --filter-min-enable 100
  "$MMSEQS" result2msa "${INDB}" "${DB3}.idx" "${BASE}/${INBASE}_res_env_ext_filter" $DB3ALI --msa-format-mode 6 --db-load-mode 2 ${FILTER_PARAM}
  "$MMSEQS" convertalis "${BASE}/${INBASE}_prof_res_env_ext" $DB3 "${BASE}/${INBASE}_res_env_ext_filter" "${OUTDB}.m8" --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar --db-load-mode 2
  "$MMSEQS" mergedbs $DB3ALI $FINALALI $MERGEALI $DB3ALI
  "$MMSEQS" rmdb "${BASE}/${INBASE}_prof_res_env_ext"
  "$MMSEQS" rmdb "${BASE}/${INBASE}_prof_res_env_ext_filter"
fi

rm -r "${LTMPDIR}"
>&2 echo $INFILE
