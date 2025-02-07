#!/bin/bash
# This script orchestrates the pipeline by submitting jobs for each FASTQ pair.
# Usage: ./bash_main_pipeline.sh INPUT_DIR OUTPUT_DIR SRC_DIR
# Example: ./bash_main_pipeline.sh /path/to/input /path/to/output /path/to/src

# Read command-line arguments for directories
INPUT_DIR="$1"
OUTPUT_DIR="$2"
src="$3"

# Create the output and run_stats directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "run_stats"

# Load required modules
module load flash/1.2.11
module load seqtk/1.0
module load cutadapt/4.4

# Export the variables so they are available to main_pipeline.sh
export INPUT_DIR
export OUTPUT_DIR
export src

# Loop through each R1 FASTQ file in the input directory
for R1 in "$INPUT_DIR"/*_R1_*.fastq.gz; do
    # Derive the R2 filename by replacing _R1_ with _R2_
    R2="${R1/_R1_/_R2_}"
    
    # Extract the basename (e.g., sample name) from the R1 filename
    BASENAME=$(basename "$R1" | cut -d'.' -f1)
    SAMPLE=${BASENAME%%_*}
    
    # Submit the job to the cluster using bsub
    # Note: We pass R1, R2, BASENAME, and SAMPLE as arguments to main_pipeline.sh
    bsub -q priority \
         -e "run_stats/${SAMPLE}.error" \
         -o "run_stats/${SAMPLE}.out" \
         -M 300000MB \
         bash "$src/main_pipeline.sh" "$R1" "$R2" "$BASENAME" "$SAMPLE"

    echo "Submitted job for $SAMPLE."
done

mkdir -p run_stats
mkdir -p plots
mkdir -p tables
mkdir -p results_log2FC
