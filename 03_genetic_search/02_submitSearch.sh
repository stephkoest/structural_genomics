#!/bin/bash
#SBATCH -J search           # Job name
#SBATCH -o log/job.%j.out   # Name of stdout output file (%j expands to jobId)
#SBATCH -e log/job.%j.err   # Name of stderr output file (%j expands to jobId)
#SBATCH -n 128              # Total number of threads or total number of mpi tasks
#SBATCH --partition=fat_genoa
#SBATCH -t 1-00:00:00

INDIR=$1

module load 2022
module load parallel/20220722-GCCcore-11.3.0

###MMSEQS only works when both modules are loaded
MMSEQS=/projects/0/lwc2020006/alpaca/software/MMseqs2/build/bin/mmseqs
VMTOUCH=/projects/0/lwc2020006/alpaca/software/vmtouch/vmtouch
DATABASES=/projects/0/lwc2020006/alpaca/databases/colabfold
DB1=${DATABASES}/uniref30_2202_db
DB2=${DATABASES}/colabfold_envdb_202108_db
DB3=${DATABASES}/asgardProteomes_db

JTHREADS=16
NPJOBS=$(expr $(expr $SLURM_NTASKS / $JTHREADS))

TARLIST=$(ls ${INDIR}/*faa.tar.gz | grep -vFf <(find results/ -name "*_search.tar.gz" | rev | cut -d"/" -f1 | rev | sed 's/_search.tar.gz//'))

for file in $(echo "${TARLIST}")
do
	tar -xzvf "${file}" -C "${TMPDIR}"
done

QLIST=$(find "${TMPDIR}" -name "*.faa")
RESLIST=$(ls -d "${TMPDIR}"/*)

"${VMTOUCH}" -l -d -t ${DB1}.idx
"${VMTOUCH}" -l -d -t ${DB2}.idx
"${VMTOUCH}" -l -d -t ${DB3}.idx

#run parallel mmseqs for unifinished proteins
parallel -j $NPJOBS "time bash bin/mmseqs_search.sh ${MMSEQS} {} ${DB1} ${DB2} ${DB3} ${JTHREADS}" ::: $QLIST
#pack up results
parallel -j $NPJOBS "tar cf - {} | pigz -9 -p ${JTHREADS} > results/{/}_search.tar.gz"  ::: $RESLIST
