#!/bin/bash
source /gpfs/work2/0/lwc2020006/alpaca/software/miniconda3/bin/activate 
conda activate /gpfs/work2/0/lwc2020006/alpaca/software/colabfold
MMSEQS=/gpfs/work2/0/lwc2020006/alpaca/software/MMseqs2/build/bin/mmseqs

INTARFILES=$@
TMPDIR=/tmp/prep
OUTDIR=$(pwd)/data
A3MD=${TMPDIR}/a3m

mkdir -p $OUTDIR $TMPDIR $TMPDIR/a3m
OUTBASE=""
for file in $(echo $INTARFILES)
do
	OUTBASE=$OUTBASE_$(basename $file .tar.gz | sed 's/_faa//' )
	pigz -p 12 -dc $file | tar -xf - -C $TMPDIR
done


for file in $(find ${TMPDIR} -name "*merged_alpaca.a3m")
do
	NAME=$(basename $(dirname $file))
	OUTFILE=${A3MD}/$(basename $file _merged_alpaca.a3m).a3m
	if [ ! -f "$OUTFILE" ]
	then
		echo working on $NAME
		#create file link for colabfold script (maybe redundant)
		LINKF=$(dirname $file)/final.a3m
		if [ -f "$OUTFILE" ]
		then
			rm ${LINKF}
		fi
		ln -s $(readlink -f $file) ${LINKF}
		colabfold_split_msas --mmseqs ${MMSEQS} $(dirname $file) .
		TMPFILE=$(head -n1 ${file} | sed 's/>//' | awk '{print $1}').a3m
		python bin/a3mRemoveSelf.py ${TMPFILE} > ${OUTFILE}
		rm ${TMPFILE}
		rm ${LINKF}
	fi
done
cd $TMPDIR

tar cf - $(basename $A3MD) | pigz -9 -p 12 > ${OUTDIR}/${OUTBASE}_search.tar.gz
rm -r $TMPDIR
