#!/bin/bash

# Set the folder path and pattern
FOLDER_PATH="results_log2FC"
PATTERN="*_combined.csv"
OUTPUT_FOLDER="results_log2FC"  # Output folder for filtered files

# Loop over each matching file and submit a job
for sample_file in "$FOLDER_PATH"/$PATTERN; do
    # Get the filename without the path
    sample_filename=$(basename "$sample_file")
    
    # Submit the job
    bsub -q normal -J "process_${sample_filename}" \
         -o "${FOLDER_PATH}/logs/${sample_filename}_%J.log" \
         "python process_results.py --sample_results_file ${sample_file} --output_folder ${OUTPUT_FOLDER}"
done
