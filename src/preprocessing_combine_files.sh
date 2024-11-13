#!/bin/bash

module load conda3/202011
folder=$1
sample_info=$2

mkdir -p "${folder}/run_stats"

bsub -q priority -e ${folder}/run_stats/preprocessing_${sample_name}.error -o ${folder}/run_stats/preprocessing_${sample_name}.out -M 800000MB -n 1 python preprocessing_combine_files.py --folder_path "$folder" --sample_info "$sample_info"
