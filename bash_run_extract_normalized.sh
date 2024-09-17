#!/bin/bash

INPUT_DIR=$1

module load conda3/202011

for fastq in "$INPUT_DIR"/*.merged.sampled*fastq; do
    SAMPLE=$(basename "$fastq" .fastq)
    bsub -q priority -e run_stats/extract_seqtik_$SAMPLE.error -o run_stats/extract_seqtik_$SAMPLE.out -M 800000MB -n 1 python extract_barcodes.py $fastq 
done