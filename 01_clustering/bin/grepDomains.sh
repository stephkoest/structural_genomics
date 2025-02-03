#!/bin/bash
cog=$1
COGASS=$2
OUTFASTA=$3
OUTD=$4

echo working $cog
grep --binary-files=text $cog $COGASS | awk '{print $2}' | sort -u > members_${cog}.tmp
grep --binary-files=text $cog $COGASS | awk '{print $2","$10"-"$11}' | sort -u > clipping_${cog}.tmp
grep -A1 -Ff members_${cog}.tmp $OUTFASTA \
    | grep -v "^--$" \
    | python bin/extractMultiDomain.py clipping_${cog}.tmp > ${OUTD}/${cog}.faa

rm members_${cog}.tmp clipping_${cog}.tmp