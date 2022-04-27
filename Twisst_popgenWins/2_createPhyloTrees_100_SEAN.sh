#!/bin/bash
#
#-------------------------------------------------------------
#example script for running a single-CPU serial job via SLURM
#-------------------------------------------------------------
#
#SBATCH --job-name=CreatePhyloTree
#SBATCH --output=out_CreatePhyloTree.out
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
export OMP_NUM_THREADS=16
#
#load the respective software module you intend to use
module load python/3.7
module load phyml/3.1

#
#
#run the respective binary through SLURM's srun
#To get neighbour joining trees for snp windows, I Martin S. wrote this script that runs Phyml for windows.
#using parallelisation, and outputs a single trees file. The script can be gound in Martin's repo genomics_general.
#https://github.com/simonhmartin/genomics_general
srun --cpu_bind=verbose python /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/phyml_sliding_windows.py -T 16 -g /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/arc_sax.geno.gz --prefix output100.phyml_bionj.w100.v2 -w 100 --windType sites --model GTR --optimise n

wait
