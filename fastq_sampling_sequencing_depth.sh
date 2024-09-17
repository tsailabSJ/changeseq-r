#!/bin/bash

dir=$1
file=$2

INPUT_COUNT=$(grep -c "^@" "$file")

# Extract the sample ID from the filename
sample_id=$(basename "$file" .merged.fastq)

# Use seqtk to randomly sample the file with different fractions and output to different files
for fraction in 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
do
    # Set the output file name based on the fraction
    output_file="${dir}/${sample_id}.merged.sampled.${fraction}.fastq"

    # Use seqtk to sample the reads
    seqtk sample -s100 "$file" $fraction > "$output_file"

    echo "Sampled $fraction of $file and saved to $output_file"
done