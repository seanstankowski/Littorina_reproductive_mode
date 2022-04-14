#!/bin/bash

IN_DIR=$1
OUT_DIR=$2

FILES=($(ls -1 ${IN_DIR}/*.bam))


# get size of array
NUMFILE=${#FILES[@]}

mkdir -p ${OUT_DIR}/logs
mkdir -p ${OUT_DIR}/tmp

# now submit to queue
if [ $NUMFILE -ge 0 ]; then
qsub -t 1-$NUMFILE -N gatk_prepare_\$TASK_ID -o ${OUT_DIR}/logs/prepare_\$TASK_ID.log -e ${OUT_DIR}/logs/prepare_\$TASK_ID.err do_prepare.sh ${IN_DIR} ${OUT_DIR} 
fi
