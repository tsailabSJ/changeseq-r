#!/bin/bash

# bash bash_run_seqtik.sh ./240711_VH00889_169_AAFY7MYM5 ./merged_fastq_seqtik

INPUT_DIR=$1
OUTPUT_DIR=$2

mkdir -p "$OUTPUT_DIR"
mkdir -p "run_stats"

module load flash/1.2.11 

for R1 in "$INPUT_DIR"/*_R1_*.fastq.gz; do
    R2="${R1/_R1_/_R2_}"
    
    BASENAME=$(basename "$R1"| cut -d'.' -f1)
    SAMPLE=${BASENAME%%_*}
    
    bsub -q standard -e run_stats/merge_seqtik_$SAMPLE.error -o run_stats/merge_seqtik_$SAMPLE.out -M 400000MB -n 1 bash merge_fastq_seqtik.sh $R1 $R2 $OUTPUT_DIR
    
done