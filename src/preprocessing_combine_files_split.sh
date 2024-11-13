#!/bin/bash

module load conda3/202011

# Define the parent directory that contains the folders
# parent_directory="results"
parent_directory=$1

# Loop through all the folders in the parent directory
for folder in "$parent_directory"/*/
do
  # Remove the trailing slash from the folder path
  folder_path="${folder%/}"
  echo $folder_path
  
  sample_name=$(basename "$folder_path")
  # Run the Python script and pass the folder path as an argument
  bsub -q priority -e run_stats/preprocessing_${sample_name}.error -o run_stats/preprocessing_${sample_name}.out -M 800000MB -n 1 python preprocessing_combine_files.py --folder_path "$folder_path"

done
