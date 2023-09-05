#!/bin/bash -ve

#$ -P littorina
#$ -q littorina.q
#$ -pe smp 8
# request memory for job (default 6G, max 72G)
#$ -l mem=16G
#$ -l rmem=16G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=336:00:00


module load apps/java

/data/bo1srs/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 8 -m GTRGAMMA -p 1689131 -s /fastdata/bo1srs/final_working_VCFs/variants_only_108_VCF_MAC2_minQ30_maxmiss1.0.recode.min4.phy -w  /fastdata/bo1srs/final_working_VCFs/RAXML/ -n variants_only_108_VCF_MAC2_minQ30_maxmiss1.0_test
