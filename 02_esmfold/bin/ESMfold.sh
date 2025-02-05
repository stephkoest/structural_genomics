#!/bin/bash
#SBATCH -J prediction       # Job name
#SBATCH -o log/job.%j.out   # Name of stdout output file (%j expands to jobId)
#SBATCH -e log/job.%j.err   # Name of stderr output file (%j expands to jobId)
#SBATCH -n 72             # Total number of threads or total number of mpi tasks
#SBATCH --partition=gpu
#SBATCH -t 5-00:00:00
#SBATCH --gpus-per-node=4

module load 2022
module load Python/3.10.4-GCCcore-11.3.0
INFASTA=$1
OUTPDB=$2

source /gpfs/work2/0/lwc2020006/alpaca/software/ESMfold/bin/activate

mkdir -p $OUTPDB
python bin/esmfold_inference.py -i $INFASTA -o $OUTPDB --num-recycles 12
