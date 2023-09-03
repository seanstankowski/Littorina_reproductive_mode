#!/bin/bash

BAM_DIR=$1
FASTA=$2
NLINES=$3

## split the fasta fai file into 1000 pieces
split ${FASTA}.fai -l ${NLINES} contigs-
ls -1 contigs-* > contigs.list


for c in contigs-*
do

tmp=$(echo ${c} | sed 's/contigs-//')
mkdir -p $tmp
mkdir -p $tmp/logs

BED_FILE=`pwd`/$tmp/$tmp.bed
NEW_FASTA=`pwd`/$tmp/$tmp.fasta

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $c > ${BED_FILE}

module load apps/java
module load apps/SAMtools
PATH=$PATH:/data/bo1srs/gatk-4.0.7.0
export PATH
# May need to edit to be path to the picard.jar that you have access to
PICARD=/data/bo1srs/picard/build/libs/picard.jar

## remove any fasta and dictionary files and make a fresh copy
rm -r ${NEW_FASTA} || true
rm -r `pwd`/$tmp/$tmp.fasta.fai || true
rm -r `pwd`/$tmp/$tmp.dict || true


# create a fasta file for the contigs if it does not already exist
# make list of contig names
cat ${BED_FILE} | cut -f1 > `pwd`/${tmp}/contigs.txt
## subset original according to list of contigs
xargs samtools faidx ${FASTA} < `pwd`/${tmp}/contigs.txt > ${NEW_FASTA}
## create indices required for GATK and Picard
java -jar ${PICARD} CreateSequenceDictionary R=${NEW_FASTA} O=`pwd`/${tmp}/${tmp}.dict TMP_DIR=`pwd`/tmp
samtools faidx ${NEW_FASTA} 

FILES=($(ls -1 ${BAM_DIR}/*.bam))


# get size of array
NUMFILE=${#FILES[@]}


# now submit to queue
if [ $NUMFILE -ge 0 ]; then
qsub -t 1-${NUMFILE} -N gatk_contig_${tmp} -o `pwd`/${tmp}/logs/\$TASK_ID.log -e `pwd`/${tmp}/logs/\$TASK_ID.err do_gatk.sh ${BAM_DIR} `pwd`/${tmp} ${NEW_FASTA} ${BED_FILE}

# line below can be used to submit 5 contigs as a test
#qsub -t 1-5 -N gatk_contig_${tmp} -o `pwd`/${tmp}/logs/\$TASK_ID.log -e `pwd`/${tmp}/logs/\$TASK_ID.err do_gatk.sh ${BAM_DIR} `pwd`/${tmp} ${NEW_FASTA} ${BED_FILE}
fi

rm $c
done
