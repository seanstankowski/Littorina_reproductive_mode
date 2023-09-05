The scripts in this directory are used to prepare bams and conduct variant calling of sequence data for a fragmented reference genome. Calls are conducted for subsets of contigs, and these these are concatenated together at the end to produce a complete VCF. 

Running the pipeline

## 1. prepare bams:
The first step removes PCR duplicates from bars and adds read group info prior to calling.
	### relevant scripts
	setup.prepare
	do.prepare
	
./setup_prepare.sh path_to_bam_dir /path_to_analysis_ready_bams /path_to_genome_fasta


## 2. produce gvcf for each individual for n sets of n contigs:
This analysis produces a gvcf for each individual, for each set of contigs. We used batches of 1000 contigs to make this runnable on our hps.
	### relevant scripts
	setup.gatk.sh
	do.gatk.sh

./setup_gatk.sh /path_to_prepared_bams /path_to_genome number_of_contigs_per_job


## 3. combine calls into a vcf
This performs joint calling across all individuals for batch of contigs. 
	setup_combine.sh
	do_combine.sh

./setup_combine.sh /path_to_prepared_bams /path_to_contig_list/filename 

## 4. Put together all combined VCFs to produce a full VCF. I did this using VCFtools concat.

