#!/bin/bash
#SBATCH -J search           # Job name
#SBATCH -o log/job.%j.out   # Name of stdout output file (%j expands to jobId)
#SBATCH -e log/job.%j.err   # Name of stderr output file (%j expands to jobId)
#SBATCH -n 64              # Total number of threads or total number of mpi tasks
#SBATCH --partition=thin
#SBATCH -t 2-00:00:00

module load 2021

SEQDB="${1}"
DB="${2}"

source '/projects/0/lwc2020006/alpaca/software/miniconda3/bin/activate'

conda activate "mmseqs2_v14"

if [[ -z "${SLURM_NTASKS}" ]]
then
	SLURM_NTASKS=10
fi

MMSEQS='/gpfs/work2/0/lwc2020006/alpaca/software/miniconda3/envs/mmseqs2_v14/bin/mmseqs'

RESDB=results/$(basename ${SEQDB})_$(basename ${DB})

"${MMSEQS}" search ${SEQDB} ${DB} $RESDB $TMPDIR -e 0.001 -s 9 --remove-tmp-files 1 --threads $SLURM_NTASKS --split-memory-limit $(expr $(expr ${SLURM_NTASKS} \* 2) - 16 )"G"

echo "done search ${SEQDB} ${DB}"
"${MMSEQS}" convertalis ${SEQDB} ${DB} $RESDB ${RESDB}.m8 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits
