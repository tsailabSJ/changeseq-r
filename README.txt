CHANGE-seq randomized pipeline

smb://jude.stjude.org/hematology/common/sequencing/240801_VH00889_179_AAFT5F7M5/Data/Intensities/BaseCalls

# 0.25, 0.50, 1 data:
ln -s ~/sequencing/240801_VH00889_179_AAFT5F7M5/Data/Intensities/BaseCalls ./240801_VH00889_179_AAFT5F7M5

# 5pg data:
ln -s ~/sequencing/240903_VH00889_189_AAFLNNJM5/Data/Intensities/BaseCalls ./240903_VH00889_189_AAFLNNJM5


Step 1a: Merge the fastq files
Run time = approx. 1.5 hr

bash bash_run_seqtik.sh ./240801_VH00889_179_AAFT5F7M5 ./merged_fastq_seqtik
bash bash_run_seqtik.sh ./240903_VH00889_189_AAFLNNJM5 ./merged_fastq_seqtik

This loops through each of the fastq R1/R2 files and runs merge_fastq_seqtik.sh
Since this is CHANGE-seq (circularized DNA), the reads will map outwards. 
To properly process this, we will reverse complement R1 and concat with seqtik. 
The output files will be in merged_fastq_seqtik. 




Step 1b: Sample the fastq files

bash bash_sequencing_depth.sh 

This file loops through every *.merged.fastq files in a folder and runs fastq_sampling.sequencing_depth.sh which samples the files at a 1.0, 0.9, 0.8, 0.7, ..., 0.1 ratio and exports with the same file name plus the ratio. 




Step 2: Extract barcodes

bash bash_run_extract_normalized.sh ./merged_fastq_seqtik

Run time = approx. 2 hr

This runs extract_barcodes.py for each *.merged.fastq file in the directory.
It requires a template_sequences.csv with the following columns: Sample (ex. AF211, AF212, etc.) and Template (ex. ACGCGAGCTGCATGTGTCAGANNWNNWNNWNNWNNWGGGGCCACTAGGGACAGGATTGGCTTCTTCAAGTCCGCCATGCCCGA) 

In the extract_barcodes.py, you may need to update the positions:

# Positions of the different parts within the template
universal_primer_start = 0
universal_primer_end = 21
barcode_start = 21
barcode_end = 36
target_start = 36
target_end = 59
flanking_start = 59
flanking_end = 83

For a read to pass alignmenet, the 3' and 5' flanking regions (10 bp) must perfectly align 90% (error rate allowed: 1 out of 10). 
Then the 15 bp barcode followed by the 23 bp target region is extracted.
The output will be in results: *_results.csv and *_poor_alignment.csv





Step 3: Preprocessing the extracted barcodes. Combine to one file. 

Do each sample separately:
bash preprocessing_combine_files_split.sh

# bsub -q priority -e run_stats/preprocessing.error -o run_stats/preprocessing.out -M 800000MB -n 1 python preprocessing_combine_files.py


Run .ipynb to analysis the number of unique barcode+target reads. 