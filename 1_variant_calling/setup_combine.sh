#!/bin/bash

BAM_DIR=$1
CONTIGS_LIST=$2

NCONTIGS=$(wc -l ${CONTIGS_LIST} | cut -d " " -f 1)

mkdir -p `pwd`/logs

# now submit to queue
if [ $NCONTIGS -ge 0 ]; then
qsub -t 1-$NCONTIGS -N combine_vcf_\$TASK_ID -o `pwd`/logs/combine_vcf_\$TASK_ID.log -e `pwd`/logs/combine_vcf_\$TASK_ID.err do_combine_vcf.sh ${BAM_DIR} ${CONTIGS_LIST}
fi
