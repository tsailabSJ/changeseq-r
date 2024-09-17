#!/bin/bash

R1=$1
R2=$2
OUTPUT_DIR=$3

BASENAME_R1=$(basename "$R1"| cut -d'.' -f1)
BASENAME_R2=$(basename "$R2"| cut -d'.' -f1)
SAMPLE=${BASENAME_R1%%_*}


# reverse complement -r -p
seqkit seq -r -p $R1 > $OUTPUT_DIR/$BASENAME_R1.rc
seqkit seq $R2 > $OUTPUT_DIR/$BASENAME_R2.rc

echo "$OUTPUT_DIR/$BASENAME_R1.rc and $OUTPUT_DIR/$BASENAME_R2.rc done."


cd $OUTPUT_DIR

# merge reads, need a lot of memory, 70M reads needs 100G
seqkit concat -j 2 $BASENAME_R1.rc $BASENAME_R2.rc > $SAMPLE.merged.fastq

echo "$BASENAME_R1.rc $BASENAME_R2.rc > $SAMPLE.merged.fastq done."

module load fastqc
fastqc $SAMPLE.merged.fastq

echo "$SAMPLE.merged.fastq quality check done."