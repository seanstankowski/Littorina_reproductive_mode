#!/bin/bash
#
#-------------------------------------------------------------
#example script for running a single-CPU serial job via SLURM
#-------------------------------------------------------------
#
#SBATCH --job-name=popGenwins
#SBATCH --output=popgenWins.log
#
#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=240:00:00
#
#Define the amount of RAM used by your job in GigaBytes
#SBATCH --mem=4G
#
#Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=sean.stankowski@ist.ac.at
#SBATCH --mail-type=ALL
#
#Pick whether you prefer requeue or not. If you use the --requeue
#option, the requeued job script will start from the beginning, 
#potentially overwriting your previous progress, so be careful.
#For some people the --requeue option might be desired if their
#application will continue from the last state.
#Do not requeue the job in the case it fails.
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#for single-CPU jobs make sure that they use a single thread
export OMP_NUM_THREADS=5
#
#load the respective software module you intend to use
module load python/3
 
srun --cpu_bind=verbose python3 /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/genomics_general-master/popgenWindows.py --windType predefined --windCoords wincoords.bed -p Egglayer -p Brooder -g arc_sax.geno.gz -o popgenOut.csv.gz -f phased -T 5 --popsFile pops.txt

wait
