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

## Usage

1. Create a symbolic link to your sequencing data directory:
   ```bash
   ln -s /path/to/sequencing/241029_VH00889_211_AACLKYWHV/Data/Intensities/BaseCalls ./241029_VH00889_211_AACLKYWHV
   ```

2. Run the main pipeline script to extract randomized targets and unique barcodes from `.fastq` files:
   ```bash
   ./bash_main_pipeline.sh INPUT_DIR OUTPUT_DIR SRC_DIR
   ```
   
   Replace `INPUT_DIR`, `OUTPUT_DIR`, and `SRC_DIR` with your actual paths.

3. Create necessary directories:
   ```bash
   mkdir -p run_stats
   mkdir -p plots
   mkdir -p tables
   mkdir -p results_log2FC
   ```

4. Once the pipeline finishes (approx. runtime: ~30 minutes), combine results into one file:
   ```bash
   bsub -q priority -e run_stats/combine.error -o run_stats/combine.out -M 800000MB -n 1 \
     python src/preprocessing_combine_files.py --folder_path results_cutadaptv2 --sample_info sample_info.csv
   ```

5. For additional statistics, load conda3 and run:
   ```bash
   module load conda3/202011

   bsub -q priority -e run_stats/unique_barcodes.error -o run_stats/unique_barcodes.out -M 800000MB \
     python plots_preprocessing_unique_barcodes.py -f results_cutadaptv2/all_results_combined_cleaned_21_25.csv

   bsub -q priority -e run_stats/collision.error -o run_stats/collision.out -M 800000MB \
     python plots_preprocessing_and_collision.py
   ```

## Collision Filtering

Filter out collisions. A collision is defined as two randomized target with the same barcode. To avoid collision, we required that the Cas9 end-repaired target is within 2 edit distances of the uncleaved control sequence.

```bash
bsub -q priority -e run_stats/results_collision_filtering.error -o run_stats/results_collision_filtering.out -M 800000MB \
  python src/removing_collision_preprocessing_pybtree.py
```

## Calculating Relative Activity

Relative activity is calculated as `(fold change off-target / fold change on-target) * 100`, where:

- Fold change = (Cas9 counts) / (Control counts).
- On-target activity is used as the normalization reference.

Before running, ensure the input folder name is correctly set and samples are properly paired (control vs. treatment).

```bash
bsub -q priority -e run_stats/results_log2FC.error -o run_stats/results_log2FC.out -M 800000MB \
  python src/results_log2FC.py
```
