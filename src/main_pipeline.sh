#!/bin/bash
# This script runs the pipeline steps for a given pair of FASTQ files.
# It relies on INPUT_DIR, OUTPUT_DIR, and src variables exported by bash_main_pipeline.sh.
# Arguments:
#   $1 = R1 filename
#   $2 = R2 filename
#   $3 = BASENAME (derived from input file name)
#   $4 = SAMPLE (derived from input file name)

# Parse arguments
R1="$1"
R2="$2"
BASENAME="$3"
SAMPLE="$4"

# INPUT_DIR, OUTPUT_DIR, and src are expected to be exported by bash_main_pipeline.sh
# so we can use them directly here.

# Step 1: Merge FASTQ files using seqtk-based script
# The merge_fastq_seqtik.sh script likely merges R1 and R2 into a single merged FASTQ
bash "$src/merge_fastq_seqtik.sh" "$R1" "$R2" "$OUTPUT_DIR"

# Optional Step 2: Calculate sequencing depths
# bash "$src/fastq_sampling_sequencing_depth.sh" "${OUTPUT_DIR}" "${BASENAME}.merged.fastq"

# Step 3: Extract barcodes or other sequences of interest
# extract_barcodes.py takes a merged FASTQ and finds sequences flanked by adapters (-aL and -aR)
python "$src/extract_barcodes.py" \
    -f "${OUTPUT_DIR}/${SAMPLE}.merged.fastq" \
    -aL ATGTGTCAGA \
    -aR CTTCTTCAAG