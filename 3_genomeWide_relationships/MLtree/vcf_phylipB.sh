#!/bin/bash -ve

#$ -P littorina
#$ -q littorina.q
#$ -pe smp 16
# request memory for job (default 6G, max 72G)
#$ -l mem=10G
#$ -l rmem=10G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=336:00:00

module load apps/python/anaconda3-4.2.0

python3 /fastdata/bo1srs/final_working_VCFs/VCF2phy.py -vcf /fastdata/bo1srs/final_working_VCFs/variants_only_108_VCF_MAC2_minQ30_maxmiss1.0.recode.vcf -o /fastdata/bo1srs/final_working_VCFs/test
