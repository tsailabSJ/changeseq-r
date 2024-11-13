#!/bin/bash
# Combining all the steps into one pipeline!

# ln -s /research/dept/hem/common/sequencing/241029_VH00889_211_AACLKYWHV/Data/Intensities/BaseCalls ./241029_VH00889_211_AACLKYWHV

INPUT_DIR=/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/users/Jackie/CHANGE-seq/2024-10-31-CHANGE-seq-randomized-Part1-1mMMg/241029_VH00889_211_AACLKYWHV/
OUTPUT_DIR=results_new_new
src=/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/users/Jackie/CHANGE-seq/2024-10-31-CHANGE-seq-randomized-Part1-1mMMg-cutadapt/src

mkdir -p "${OUTPUT_DIR}"
mkdir -p "run_stats"

module load flash/1.2.11
module load seqtk/1.0
module load cutadapt/4.4

for R1 in $INPUT_DIR/*_R1_*.fastq.gz; do
    R2="${R1/_R1_/_R2_}"
    BASENAME=$(basename "$R1" | cut -d'.' -f1)
    SAMPLE=${BASENAME%%_*}
    # Submit the job using bsub
    bsub -q standard -e "run_stats/${SAMPLE}.error" -o "run_stats/${SAMPLE}.out" -M 300000MB bash $src/main_pipeline.sh $R1 $R2 $BASENAME $SAMPLE

    echo "Submitted job for $SAMPLE."
done
