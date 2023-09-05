#### Littorina haplotype blocks ####

## Load modules
module load bcftools bedops R
source ~/.bashrc
export PYTHONPATH=$HOME/.local/bin

## Working Directory
cd ~/Littorina/HaploBlocks

## Set main variables
WORKDIR=`pwd`


#-----

## Manipulate samples from .csv to usable format for ARGweaver, each sampleID needs
## to be changed to 2 haplotypes with suffixes _1 and _2
tail +2 samples_Haplo.csv | cut -f1 -d, > samples_Haplo.txt

# make sample list for egg-layers(arcana + compressa) and brooders(saxatilis)
tail +2 samples_Haplo.csv | grep -v saxatilis | cut -f1 -d,  > egg.txt # all L.arcana and compressa
tail +2 samples_Haplo.csv | grep saxatilis | cut -f1 -d,  > brood.txt # all L.saxatilis
# fix one sample name manually W_arc_04_Lamerged_sorted.bam


## Make symbolic link to the vcf file. 
mkdir VCFs && cd "$_"
ln -s  $WORKDIR/../VCFs/PHASED_variants_only_Maxmiss0.1_softDP10_minDP1500_all_sites_108_littorina.vcf.gz.vcf.gz full.vcf.gz
tabix full.vcf.gz

ln -s ~/Littorina/SWplot/VCFs_old/snail_VCFs/PHASED_all_sites_108_littorina.vcf.gz all_sites.vcf.gz
tabix all_sites.vcf.gz

## Check all is working
bcftools view -H full.vcf.gz  | head -1 # first line of the VCF file
bcftools query -l full.vcf.gz # list all sample names  

bcftools view -H full.vcf.gz | wc -l # 18539995 SNPs (variant sites)
bcftools view -H all_sites.vcf.gz | wc -l #

## No. of variants (remember to grep whole words) - each contig is indexed at 1.
bcftools view -H full.vcf.gz| grep -w -c Contig1808 # 1761 SNPs 
bcftools view -H -r Contig1808 full.vcf.gz | wc -l # 1761 SNPs  
bcftools view -H -r Contig1808 full.vcf.gz | tail -1 # last SNP at pos - 85918

bcftools view -H full.vcf.gz| grep -w -c Contig3201 # 1865 SNPs 
bcftools view -H -r Contig3201 full.vcf.gz | wc -l # 1865 SNPs
bcftools view -H -r Contig3201 full.vcf.gz | tail -1 # last SNP at pos - 101393

## Subset individuals and region from the full vcf
# check all individuals names
bcftools query -l full.vcf.gz
#fix this name: W_arc_04_Lamerged_sorted.bam_B


#+++++
######
###### Contig1808 ######

## Set variables
chrom=Contig1808
chromStart=1
chromEnd=30000

## Subset indiviudals and region
bcftools view -Oz -r ${chrom}:${chromStart}-${chromEnd} -S $WORKDIR/sample/samples_Haplo.txt  full.vcf.gz > ${chrom}_Haplo.vcf.gz
tabix ${chrom}_Haplo.vcf.gz

bcftools view -H ${chrom}_Haplo.vcf.gz | wc -l # 577 SNPs
bcftools view -H ${chrom}_Haplo.vcf.gz | tail -1 # last SNP at pos - 29957

## Change to sites format
bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' ${chrom}_Haplo.vcf.gz > ${chrom}_Haplo.gt
# USAGE: Rscript vcf2sites.R <sample.txt> <input.gt> <chr> <begin> <end> #
Rscript $WORKDIR/scripts/vcf2sites.R $WORKDIR/sample/samples_Haplo.txt ${chrom}_Haplo.gt $chrom $chromStart $chromEnd
cat <(head -2 ${chrom}_Haplo_RefAlt.sites) <(cut -f1,2 <(tail +3 ${chrom}_Haplo_RefAlt.sites)) > ${chrom}.sites

## Run ARGweaver unmasked
mkdir ${chrom}_run1 && cd $_

~/.local/bin/arg-sample -s ../VCFs/Contig1808.sites \
	-N 2e5 \
	-m 1.2e-8 \
	-r 1.2e-8 \
	--ntimes 30 \
	--maxtime 4e6 \
	--iters 6000 \
	--delta 0.01 \
	--output sample/${chrom}

sbatch ../scripts/slurm/job-runArgweaver.slurm ../VCFs/Contig1808.sites sample/Contig1808

## make mask file
# subset Contig from all_sites.vcf
bcftools view -Oz -r ${chrom}:${chromStart}-${chromEnd} -S $WORKDIR/sample/samples_Haplo.txt  all_sites.vcf.gz > ${chrom}_all_sites.vcf.gz
tabix ${chrom}_all_sites.vcf.gz

bcftools view -H ${chrom}_all_sites.vcf.gz | wc -l # 13948 sites (variants + invariants)

# mask file
bcftools view -H ${chrom}_all_sites.vcf.gz | awk '$0 !~ /^#/ {print $1,$2-1,$2}' | bedops -d <(echo -e "$chrom\t$(($chromStart-1))\t$chromEnd") - | bgzip > ${chrom}_mask.bed.gz

## Run ARGweaver masked
mkdir ${chrom}_masked_run1 && cd $_

~/.local/bin/arg-sample -s $1 \
	--maskmap $2
	-N 2e5 \
	-m 1.2e-8 \
	-r 1.2e-8 \
	--ntimes 30 \
	--maxtime 4e6 \
	--iters 6000 \
	--delta 0.01 \
	--output $3

sbatch ../scripts/slurm/job-runArgweaver_masked.slurm ../VCFs/Contig1808.sites ../VCFs/Contig1808_mask.bed.gz sample/Contig1808


## Post-process
smc2bed-all ./sample/${chrom}

burnin=5000 #decide burnin by looking at the MCMC traces

## Estimate TMRCA for all individuals, only egg-layers and brooders
for pop in egg brood samples_Haplo; do
    echo Summarizing $pop
    arg-summarize --arg-file ./sample/$chrom.bed.gz \
			 --subset-inds ../sample/$pop.txt \
			 --log-file ./sample/$chrom.log \
			 --tmrca \
			 --burnin $burnin \
			 --mean --stdev --quantile 0.05,0.5,0.95 |\
		bgzip > $chrom.tmrca.$pop.bed.gz
		echo " -> Done"
done

# Do the same for the masked by changing the directory. 


chrom=Contig1808
iter=10000
gunzip ${chrom}.${iter}.smc.gz
tail +3 $chrom.${iter}.smc | grep TREE > $chrom.${iter}_trees.txt
tail +3 $chrom.${iter}.smc | grep SPR > $chrom.${iter}_SPR.txt

## Estimate tmrca for all populations only for iteration 10000
for pop in egg brood samples_Haplo; do
    echo Summarizing $pop
    arg-summarize --arg-file ./sample/$chrom.bed.gz \
			 --subset-inds ../sample/$pop.txt \
			 --log-file ./sample/$chrom.log \
			 --tmrca \
			 --sample $iter \
			 --mean --quantile 0.05,0.5,0.95 |\
		bgzip > $chrom.${iter}.tmrca.$pop.bed.gz
		echo " -> Done"
done



## Get names from the argweaevr IDs
head -1 Contig1808.10000.smc | tr "\t" "\n" > ../../sample/argID_hap.txt

###### contig3201 unmasked

# set variable and repeat the same codes above.
chrom=Contig3201
chromStart=45000
chromEnd=75000

## Subset indiviudals and region
bcftools view -Oz -r ${chrom}:${chromStart}-${chromEnd} -S $WORKDIR/sample/samples_Haplo.txt  full.vcf.gz > ${chrom}_Haplo.vcf.gz
tabix ${chrom}_Haplo.vcf.gz

bcftools view -H ${chrom}_Haplo.vcf.gz | wc -l # 847 SNPs
bcftools view -H ${chrom}_Haplo.vcf.gz | tail -1 # last SNP at pos - 74740

## Change to sites format
bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' ${chrom}_Haplo.vcf.gz > ${chrom}_Haplo.gt
# USAGE: Rscript vcf2sites.R <sample.txt> <input.gt> <chr> <begin> <end> #
Rscript $WORKDIR/scripts/vcf2sites.R $WORKDIR/sample/samples_Haplo.txt ${chrom}_Haplo.gt $chrom $chromStart $chromEnd
cat <(head -2 ${chrom}_Haplo_RefAlt.sites) <(cut -f1,2 <(tail +3 ${chrom}_Haplo_RefAlt.sites)) > ${chrom}.sites

## Run ARGweaver without mask
mkdir ${chrom}_run1 && cd $_

sbatch ../scripts/slurm/job-runArgweaver.slurm ../VCFs/Contig3201.sites sample/Contig3201

## make mask file
# subset Contig from all_sites.vcf
bcftools view -Oz -r ${chrom}:${chromStart}-${chromEnd} -S $WORKDIR/sample/samples_Haplo.txt  all_sites.vcf.gz > ${chrom}_all_sites.vcf.gz
tabix ${chrom}_all_sites.vcf.gz

bcftools view -H ${chrom}_all_sites.vcf.gz | wc -l # 17021 sites (variants + invariants)

# mask file
bcftools view -H ${chrom}_all_sites.vcf.gz | awk '$0 !~ /^#/ {print $1,$2-1,$2}' | bedops -d <(echo -e "$chrom\t$(($chromStart-1))\t$chromEnd") - | bgzip > ${chrom}_mask.bed.gz

## Run ARGweaver masked
mkdir ${chrom}_masked_run1 && cd $_

~/.local/bin/arg-sample -s $1 \
	--maskmap $2
	-N 2e5 \
	-m 1.2e-8 \
	-r 1.2e-8 \
	--ntimes 30 \
	--maxtime 4e6 \
	--iters 6000 \
	--delta 0.01 \
	--output $3

sbatch ../scripts/slurm/job-runArgweaver_masked.slurm ../VCFs/Contig3201.sites ../VCFs/Contig3201_mask.bed.gz sample/Contig3201

## Resume ARGweaver till 10000
~/.local/bin/arg-sample -s $1 \
	--maskmap $2
	-N 2e5 \
	-m 1.2e-8 \
	-r 1.2e-8 \
	--ntimes 30 \
	--maxtime 4e6 \
	--iters $3 \
	--delta 0.01 \
	--output $4 \
	--resume
#put iters as an input in the sbatch command


## Post-process
smc2bed-all ./sample/${chrom}

burnin=5000 #decide burnin by looking at the MCMC traces

## Estimate TMRCA for all individuals, only egg-layers and brooders
for pop in egg brood samples_Haplo; do
    echo Summarizing $pop
    arg-summarize --arg-file ./sample/$chrom.bed.gz \
			 --subset-inds ../sample/$pop.txt \
			 --log-file ./sample/$chrom.log \
			 --tmrca \
			 --burnin $burnin \
			 --mean --stdev --quantile 0.05,0.5,0.95 |\
		bgzip > $chrom.tmrca.$pop.bed.gz
		echo " -> Done"
done

# Do the same for the masked by changing the directory. 

prefix=Contig3201
iter=10000
gunzip ${prefix}.${iter}.smc.gz
tail +3 $prefix.${iter}.smc | grep TREE > $prefix.${iter}_trees.txt
tail +3 $prefix.${iter}.smc | grep SPR > $prefix.${iter}_SPR.txt

## Get names from the argweaevr IDs
head -1 ${prefix}.10000.smc | tr "\t" "\n" > ../../sample/${prefix}_argID_hap.txt

