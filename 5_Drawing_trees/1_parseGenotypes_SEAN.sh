#!/bin/bash
#
#-------------------------------------------------------------
#example script for running a single-CPU serial job via SLURM
#-------------------------------------------------------------
#
#SBATCH --job-name=parseGenotypeFile
#SBATCH --output=out_parseGenotypeFile.log
#
#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=04:10:00
#
#Define the amount of RAM used by your job in GigaBytes
#SBATCH --mem=4G
#
#Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=dgarciac@ist.ac.at
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
export OMP_NUM_THREADS=1
#
#load the respective software module you intend to use
module load python/3
#
#
#run the respective binary through SLURM's srun
#The script to generate the trees takes a simple genotype format as input.
#The script parseVCF.py removes unnecesary informatin keeping only the contig, position, and genotypes of the samples.
srun --cpu_bind=verbose python /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/parseVCF.py -i /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/full_Q40_maxmiss1_mac1_PHASED.vcf.gz --skipIndels | gzip > arc_sax.geno.gz

