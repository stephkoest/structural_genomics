#!/bin/bash
TEMP=/tmp
for file in results/Asgard_DB_remaining*_c*_pred.tar.gz
do
        echo $file
        tar -xzf ${file} -C $TEMP  && cp ${TEMP}/predictions/log.txt results/$(basename $file .tar.gz).log
        mv ${TEMP}/predictions/*rank_1*.pdb results/Asgard_COLAB
done
find results/Asgard_COLAB/ -name "*.pdb" -exec awk -v prot=$(basename {} .pdb) 'NR>1{sum += $11}END{print prot"\t"sum/NR}' {} \; > results/Asgard_COLAB_pLDDT.tsv

mkdir results/Asgard_best_structures/
#Now pick colabfold predicted structures above threshold (by default winners)
ln -s $(awk '$2>=80{print "/projects/0/lwc2020006/alpaca/prediction/"$1}' results/Asgard_COLAB_pLDDT.tsv ) results/Asgard_best_structures/
