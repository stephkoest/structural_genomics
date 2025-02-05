#!/bin/bash

NAME=Asgard_DB_remaining
#input tmp dir path
TEMP=/tmp
INFA=data/Asgard_DB/Asgard_DB_remaining_representatives.faa

mkdir -p $TEMP
for prot in $(grep ">" $INFA | sed 's/>//')
do
        grep -A1 $prot $INFA > ${TEMP}/${prot}.faa
done

LFILE=length.tmp
for file in $(find $TEMP -name "*.faa")
do
        length=$(head -n2 $file | tail -n1 | wc -m)
        echo $file $length
done > $LFILE
awk '$2 <= 300 {print $1}' $LFILE > small.tmp
awk '$2 > 300 && $2 <= 500 {print $1}' $LFILE > medium.tmp
awk '$2 > 500 && $2 <= 800 {print $1}' $LFILE > large.tmp
awk '$2 > 800 {print $1}' $LFILE > verylarge.tmp
split -d -l 3000 small.tmp small_c
split -d -l 1000 medium.tmp medium_c
split -d -l 200 large.tmp large_c
split -d -l 50 verylarge.tmp verylarge_c

WORKD=$(pwd)
for file in  *_c[0-9]*
do
        echo $file
        OUTNAME=${NAME}_${file}
        mkdir -p $TEMP/${OUTNAME}
        mv $(cat $file) $TEMP/${OUTNAME}/
        cd $TEMP
        mkdir -p ${WORKD}/data/Asgard_DB/${OUTNAME}
        tar cf - ${OUTNAME} | pigz -9 -p 8 > ${WORKD}/data/Asgard_DB/${OUTNAME}/${OUTNAME}_faa.tar.gz
        cd $WORKD
done
rm *_c[0-9]* small.tmp medium.tmp large.tmp verylarge.tmp
