#!/bin/bash
#$ -l rmem=16G
#$ -m bea # send mails at beginning, end and if aborted unexpectedly
#$ -M 
IN_DIR=$1
OUT_DIR=$2


FILES=($(ls -1 ${IN_DIR}/*.bam))

FILENAME=${FILES[$SGE_TASK_ID-1]}
sample_name=$(echo $(basename $FILENAME) | sed 's/_merged_sorted.bam//')

echo The input file is ${FILENAME}
echo The sample name is ${sample_name}

module load apps/java
PATH=$PATH:/shared/bioinformatics_core1/Shared/software/bin
export PATH
PICARD=/shared/bioinformatics_core1/Shared/software/picard.jar

if [ ! -f ${OUT_DIR}/${sample_name}.done ]
then

#java '-Xmx80g' -jar $PICARD MarkDuplicates I=${FILENAME} O=${OUT_DIR}/${sample_name}.dups.bam M=${OUT_DIR}/tmp/${sample_name}.dupMetrics.txt VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true ASSUME_SORT_ORDER=coordinate TMP_DIR=${OUT_DIR}/tmp/

bammarkduplicates index=1 I=${FILENAME} O=${OUT_DIR}/${sample_name}.dups.bam

java -jar $PICARD AddOrReplaceReadGroups \
      I=${OUT_DIR}/${sample_name}.dups.bam \
      O=${OUT_DIR}/${sample_name}.RG.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
RGSM=${sample_name} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

rm ${OUT_DIR}/${sample_name}.dups.bam
rm ${OUT_DIR}/${sample_name}.dups.bam.bai
touch ${OUT_DIR}/${sample_name}.done
else echo "Sample is ready for analysis"
fi



