# CHANGE-seq-R

This open source pipeline processes CHANGE-seq-R sequencing reads, extracts the randomized target region, and calculates the relative activity of each target.

## Overview

The pipeline:

1. Uses `cutadapt` to trim flanking sequences `ATGTGTCAGA` and `CTTCTTCAAG` from the sequencing reads.
2. Extracts a 15 bp barcode and a target region from each read.
3. Ensures the barcode matches an `NNWNNWNNWNNWNNW` pattern (15 bp total, where `W` = A/T and `N` = A/C/G/T).
4. Retains targets between 21 and 25 bp long to allow for up to 2 edit distances in downstream analysis.
5. Calculates relative activity by comparing fold changes (Cas9 vs. control and then normalizing to on-target activity).

## Adapter Trimming Details

- 10 bp adapters are used.
- Allowed error rate: 0.1.
- After trimming:
  - The first 15 bp are extracted as a barcode.
  - The remainder is considered the target region.
  
## Requirements

- **SeqKit v2.5.1**  
- **Cutadapt v4.4**

Ensure these tools are installed and accessible in your environment.

## Instructions for demo with test files
Test fastq.gz and corresponding sample.txt files are provided for testing in `test_data`. These are sampled CTLA4 and CTLA4_control R1 and R2 reads. 

1. **Barcode and radnomized target extraction:** Run the main pipeline script to merge R1 and R2 fastq.gz files using seqkit and extract the randomized targets and unique barcodes. 
   ```bash
   ./bash_main_pipeline.sh test_data results src
   ```
   The results are saved in a new folder called `results`. The results include `*__results.csv` files with two columns: barcode and randomized target for each paired samples. Reads that failed filter (i.e. the barcode that does not match the pattern or the randomized target is too long or too short) are stored in `*__failed_reads.csv`. Quality check in the form of `*fastqc.html` are also provided.

2. **Combine results:** After extracting the barcodes and randomized targets from the CHANGE-seq-R data, compare the Cas9-treated samples with their corresponding controls. Because end-repair may occur following Cas9 editing, we rely on the barcode and untreated controls to identify the original randomized target sequence. If multiple samples are present, they can be merged into a single file. 

  A `sample_info.csv` file containing sample IDs, treatment types, and additional annotations is provided to facilitate accurate comparison and integration of the samples in this test set.

   ```bash
   bsub -q priority -e run_stats/combine.error -o run_stats/combine.out -M 20000MB -n 1 \
     python src/preprocessing_combine_files.py --folder_path results --sample_info test_data/sample_info.csv
   ```

  The output will be saved in the same results folder: 
  - `all_results_combined.csv`: contains all the reads (barcode and randomized target) for each sample.
  - `all_results_combined_cleaned_21_25.csv`: contains reads that remain after filtering based on randomized target lengths

3. **Analyze results:** Read counts for each randomized target are normalized by counts per million (CPM). Then fold-change is calculated. Additonal downstream analysis can be conducted with these read counts, such as calculating the relative activity to on-target. 

  ```bash
  bsub -q priority -e run_stats/results_log2FC.error -o run_stats/results_log2FC.out -M 20000MB \
python src/results_log2FC.py --combined_results results/all_results_combined_cleaned_21_25.csv
  ```
  The output files are located in `results_log2FC/results_*_combined.csv`. This file contains the following information:
  - control_target_seq: the original unedited randomized target found in unedited controls
  - control_counts: the number of CHANGE-seq-R cleaved reads that contained this randomized target
  - Cas9_counts: the number of CHANGE-seq-R cleaved reads that contained this randomized target
  - MM: the number of mismatches between this randomized target and the sgRNA or on-target sequence
  - control_ratios: control_counts/total_control_counts
  - Cas9_ratios: Cas9_counts/total_Cas9_control_counts
  - control_norm: normalized control reads (CPM)
  - Cas9_norm: normalized Cas9 reads (CPM)
  - log2FC: log2 fold change of Cas9_norm/control_norm


## Advanced Usage

1. Create a symbolic link to your sequencing data directory:
   ```bash
   ln -s /path/to/sequencing/241029_VH00889_211_AACLKYWHV/Data/Intensities/BaseCalls ./241029_VH00889_211_AACLKYWHV
   ```

2. Run the main pipeline script to extract randomized targets and unique barcodes from `.fastq` files:
   ```bash
   ./bash_main_pipeline.sh INPUT_DIR OUTPUT_DIR SRC_DIR
   ```
   
   Replace `INPUT_DIR`, `OUTPUT_DIR`, and `SRC_DIR` with your actual paths.

3. Once the pipeline finishes (approx. runtime: ~30 minutes), combine results into one file:
   ```bash
   bsub -q priority -e run_stats/combine.error -o run_stats/combine.out -M 800000MB -n 1 \
     python src/preprocessing_combine_files.py --folder_path results_cutadaptv2 --sample_info sample_info.csv
   ```

4. For additional statistics, load conda3 and run:
   ```bash
   module load conda3/202011

   bsub -q priority -e run_stats/unique_barcodes.error -o run_stats/unique_barcodes.out -M 800000MB \
     python plots_preprocessing_unique_barcodes.py -f results_cutadaptv2/all_results_combined_cleaned_21_25.csv

   bsub -q priority -e run_stats/collision.error -o run_stats/collision.out -M 800000MB \
     python plots_preprocessing_and_collision.py
   ```

## Collision Filtering (Optional)

Filter out collisions. A collision is defined as two randomized target with the same barcode. To avoid collision, we required that the Cas9 end-repaired target is within 2 edit distances of the uncleaved control sequence.

```bash
bsub -q priority -e run_stats/results_collision_filtering.error -o run_stats/results_collision_filtering.out -M 800000MB \
  python src/removing_collision_preprocessing_pybtree.py
```

## Calculating Relative Activity (Additional data preprocessing)

Relative activity is calculated as `(fold change off-target / fold change on-target) * 100`, where:

- Fold change = (Cas9 counts) / (Control counts).
- On-target activity is used as the normalization reference.

Before running, ensure the input folder name is correctly set and samples are properly paired (control vs. treatment). This output only contains normalized read counts and log2FCs. To calculate the relative activity, you would need the on-target sequences. 

```bash
bsub -q priority -e run_stats/results_log2FC.error -o run_stats/results_log2FC.out -M 800000MB \
  python src/results_log2FC.py
```
