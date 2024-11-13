#!/bin/bash

INPUT_DIR=/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/users/Jackie/CHANGE-seq/CHANGE-seq/2024-10-31-CHANGE-seq-randomized-Part1-1mMMg/241029_VH00889_211_AACLKYWHV
OUTPUT_DIR=results_new_new
src=/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/users/Jackie/CHANGE-seq/2024-10-31-CHANGE-seq-randomized-Part1-1mMMg-cutadapt/src

R1=$1
R2=$2
BASENAME=$3
SAMPLE=$4

# # Step 1: Merge fastq seqtik
# bash "$src/merge_fastq_seqtik.sh" "$R1" "$R2" "$OUTPUT_DIR" #it should be $INPUT_DIR_merged_fastq_seqtik (and make sure there is directory)

# # Step 2: Calculate sequencing depths
# bash "$src/fastq_sampling_sequencing_depth.sh" "${INPUT_DIR}_merged_fastq_seqtik" "${BASENAME}.merged.fastq"

# Step 3: Extract sequences
python "$src/extract_barcodes.py" -f "${OUTPUT_DIR}/${SAMPLE}.merged.fastq" -aL ATGTGTCAGA -aR CTTCTTCAAG