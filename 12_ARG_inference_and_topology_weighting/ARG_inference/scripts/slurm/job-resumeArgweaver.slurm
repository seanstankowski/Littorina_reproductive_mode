#!/bin/bash

##### USAGE:sbatch job-resumeArgweaver.slurm <SITES file> <iters> <out prefix>
##### out_prefix=.../sample/...

#SBATCH --job-name=Littorina_Haplo
#SBATCH --output=Littorina.slurm-%A.out
#SBATCH --error=Littorina.slurm-%A.err

#SBATCH --partition=defaultp
#SBATCH --time=240:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=50G
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=ALL

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## load any module you need here
module load python2 samtools bedops bcftools

## run commands on SLURM's srun
srun ~/.local/bin/arg-sample -s $1 \
        -N 2e5 \
        -m 1.2e-8 \
        -r 1.2e-8 \
        --ntimes 30 \
        --maxtime 4e6 \
        --iters $2 \
        --delta 0.01 \
        --output $3 \
	--resume
