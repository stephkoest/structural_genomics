#!/bin/bash
#SBATCH -J prediction       # Job name
#SBATCH -o log/job.%j.out   # Name of stdout output file (%j expands to jobId)
#SBATCH -e log/job.%j.err   # Name of stderr output file (%j expands to jobId)
#SBATCH -n 18             # Total number of threads or total number of mpi tasks
#SBATCH --partition=gpu
#SBATCH -t 5-00:00:00
#SBATCH --gpus-per-node=1

MMSEQS=/gpfs/work2/0/lwc2020006/alpaca/software/MMseqs2/build/bin/mmseqs

source /gpfs/work2/0/lwc2020006/alpaca/software/miniconda3/bin/activate
conda activate /gpfs/work2/0/lwc2020006/alpaca/software/colabfold

INFILE=$1
tar -xzvf $INFILE -C $TMPDIR

INDIR=$(dirname $(find $TMPDIR -name "*.a3m" | head -n1))

mkdir $TMPDIR/predictions
RESD=$(pwd)/results



colabfold_batch --stop-at-score 85 --stop-at-score-below 50 $INDIR $TMPDIR/predictions
cd $TMPDIR
tar cf - predictions | pigz -9 -p 18 > ${RESD}/$(basename $INFILE .tar.gz)_pred.tar.gz
