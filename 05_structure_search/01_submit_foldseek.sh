#!/bin/bash
#SBATCH -J st_search        # Job name
#SBATCH -o log/job.%j.out   # Name of stdout output file (%j expands to jobId)
#SBATCH -e log/job.%j.err   # Name of stderr output file (%j expands to jobId)
#SBATCH -n 128               # Total number of threads or total number of mpi tasks
#SBATCH --partition=thin
#SBATCH -t 5-00:00:00

VMTOUCH=/projects/0/lwc2020006/alpaca/software/vmtouch/vmtouch
module load 2021
module load CMake/3.20.1-GCCcore-10.3.0
module load GCC/10.3.0
module load parallel/20210622-GCCcore-10.3.0

source /gpfs/work2/0/lwc2020006/alpaca/software/miniconda3/bin/activate
conda activate foldseek

#export PATH=/projects/0/lwc2020006/alpaca/software/foldseek/bin:$PATH

DB2=data/Asgard_db
MAXS=10000

mkdir -p results

for file in $(echo ../databases/foldseek/PDB_db  ../databases/foldseek/UniProt50_Alphafold_db  ../databases/foldseek/SwissProt_Alphafold_db);
do
        DB1=$file
        ALN1=results/$(basename $DB1 _db)_$(basename $DB2 _db)_aln
        ALN2=results/$(basename $DB2 _db)_$(basename $DB1 _db)_aln
        #DB1 vs DB2
        foldseek search $DB1 $DB2 $ALN1 $TMPDIR --threads $SLURM_NTASKS -a --max-seqs $MAXS
        foldseek convertalis $DB1 $DB2 $ALN1 ${ALN1}.m8 --format-mode 4 \
                --format-output query,target,theader,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits --threads $SLURM_NTASKS
        #DB2 vs DB1
        foldseek search $DB2 $DB1 $ALN2 $TMPDIR --threads $SLURM_NTASKS -a --max-seqs $MAXS
        foldseek convertalis $DB2 $DB1 $ALN2 ${ALN2}.m8 --format-mode 4 \
                --format-output query,target,theader,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,taxid,taxname,taxlineage
done
