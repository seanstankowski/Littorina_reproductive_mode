#!/bin/bash
#
#-------------------------------------------------------------
#example script for running a single-CPU serial job via SLURM
#-------------------------------------------------------------
#
#SBATCH --job-name=TwisstWeighting
#SBATCH --output=out_TwisstWeighting.log
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
module load python/3.7.9
module load phyml/3.1

#
#
#run the respective binary through SLURM's srun
#The main script, twisst.py implements the topology weighting.
#It requires ete3 and numpy libraries already installed on phyton.
#The main input is a tree file containing one or more trees in newick format. 
#Taxa (groups) must be specified in the command line, using the -g flag. This flag must be present at least four times (with three groups there is only one possible unrooted topology). 
#Alternatively (or additionally), a tab-delimited file of tip labels and their corresponding groups can be provided with the --groupsFile flag. 
#This should have tip labels in the first column and group names in the second column. 
srun --cpu_bind=verbose python /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/twisst/twisst.py -t /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/output.phyml_bionj.w100.trees.gz -w output.weights.csv.gz --outputTopos topologies.trees -g arcana -g compressa -g NorthAtlantic -g Spain --method fixed --groupsFile /nfs/scistore03/bartogrp/sstankow/TopologyWeightingFiles/TwisstScripts/groups.txt

wait
