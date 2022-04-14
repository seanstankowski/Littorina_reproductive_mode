#!/bin/bash
#$ -l rmem=32G
#$ -P littorina
#$ -m bea # send mails at beginning, end and if aborted unexpectedly
#$ -M s.stankowski@Sheffield.ac.uk # mail sent to this address

module load apps/java

BAM_DIR=$1
CONTIGS_LIST=$2

#echo "Task id is $SGE_TASK_ID"
CONTIG=$(awk "NR==$SGE_TASK_ID" $CONTIGS_LIST)

FILES=($(ls -1 ${BAM_DIR}/*.bam))
NUMFILE=${#FILES[@]}

tmp=$(echo $CONTIG | sed 's/contigs-//')

echo Directory is ${tmp}

sample_list=${tmp}/input.list

## Make sure the sample is empty
rm -r ${sample_list} || true

touch ${sample_list}

for bam in ${BAM_DIR}/*.bam
do
sample_name=$(echo $(basename $bam) | sed 's/.RG.bam//')
vcf=${tmp}/$sample_name.g.vcf.gz
echo -e "${sample_name}\t${vcf}" >> ${sample_list}
done


COMP_FILES=($(ls -1 ${tmp}/*.done))
NUMCOMP=${#COMP_FILES[@]}

echo Number of expected files is ${NUMFILE}
echo Number of completed files is ${NUMCOMP}

##The database directory needs to empty, or be created by GATK
rm -rf ${tmp}/db/ || true


if [ $NUMCOMP -eq $NUMFILE ];
then
	echo "GATK has completed for samples. Merging the vcfs..."
	module load apps/java
	PATH=$PATH:/data/bo1srs/gatk-4.2.5.0
	export PATH
	gatk --java-options "-Xmx24g -Xms24g" GenomicsDBImport --genomicsdb-workspace-path ${tmp}/db/ -L ${tmp}/${tmp}.bed --sample-name-map ${sample_list}
	gatk --java-options "-Xmx20g" GenotypeGVCFs \
	   -R ${tmp}/${tmp}.fasta \
   	   -V gendb://${tmp}/db/ \
   	   -O ${tmp}/${tmp}.combined.vcf.gz \
   	   --include-non-variant-sites \
   	   --allow-old-rms-mapping-quality-annotation-data

	status=$?
	if [ $status == 0 ];
	then
	touch ${tmp}/${tmp}.combine_vcf.done
	else
	touch ${tmp}/${tmp}.combine.vcf.fail
	fi

else
	echo "Not all GATK runs completed. Please re-run"
	for bam in ${BAM_DIR}/*.bam
	do
	echo $bams
	sample_name=$(echo $(basename $bam) | sed 's/.RG.bam//')
	echo "Expecting to find file ${tmp}/$sample_name.done"
	if [ ! -f ${tmp}/${sample_name}.done ]
	then
	echo "Analysis for ${sample_name} needs to be re-run"
	else echo "Analysis for ${sample_name} has been completed"
	fi
	done
fi

