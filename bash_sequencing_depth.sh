#!/bin/bash

module load seqtk/1.0

# Loop through each .merged.fastq file in the directory
dir='merged_fastq_seqtik'

for file in ${dir}/*.merged.fastq
do
    SAMPLE=$(basename "$file" .merged.fastq)
    bsub -q priority -e run_stats/sampling_$SAMPLE.error -o run_stats/sampling_$SAMPLE.out -M 800000MB -n 1 bash fastq_sampling_sequencing_depth.sh $dir $file
done