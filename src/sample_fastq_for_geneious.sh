#!/bin/bash

# Function to extract the first N read IDs from the CSV file, regardless of the column position
extract_read_ids() {
    csv_file=$1
    num_rows=$2
    output_file=$3

    # Determine the column position of read_id
    read_id_col=$(head -n 1 "$csv_file" | tr ',' '\n' | grep -n 'read_id' | cut -d: -f1)

    # Extract the read_ids based on the column number and store them in the output file
    awk -F, -v col="$read_id_col" 'NR>1 {print $col}' "$csv_file" | head -n "$num_rows" > "$output_file"
}

# Function to sample the FASTQ file based on read IDs
sample_fastq_by_read_id() {
    fastq_file=$1
    read_ids_file=$2
    output_fastq=$3

    # Use seqtk to sample the FASTQ file based on the list of read IDs
    seqtk subseq "$fastq_file" "$read_ids_file" > "$output_fastq"
}

# Main function
main() {
    csv_file=$1
    num_rows=$2
    fastq_file=$3

    # Create temporary files to store the read IDs and output FASTQ
    read_ids_file=$(mktemp)
    output_fastq="${csv_file%.csv}.sampled.fastq"

    # Extract read IDs
    extract_read_ids "$csv_file" "$num_rows" "$read_ids_file"

    # Sample the FASTQ file using the extracted read IDs
    sample_fastq_by_read_id "$fastq_file" "$read_ids_file" "$output_fastq"

    # Clean up the temporary read_ids file
    rm -f "$read_ids_file"

    echo "Sampled FASTQ file saved as: $output_fastq"
}

# Example usage
num_rows=5000
fastq_file="../240930_VH00889_198_AAC527NHV_merged_fastq_seqtik/AF237.merged.sampled.0.1.fastq"

csv_file="../results_sampling/AF237_0.1_poor_alignments.csv"
main "$csv_file" "$num_rows" "$fastq_file"

csv_file="../results_sampling/AF237_0.1_results.csv"
main "$csv_file" "$num_rows" "$fastq_file"




fastq_file="../240930_VH00889_198_AAC527NHV_merged_fastq_seqtik/AF249.merged.sampled.0.1.fastq"

csv_file="../results_sampling/AF249_0.1_poor_alignments.csv"
main "$csv_file" "$num_rows" "$fastq_file"

csv_file="../results_sampling/AF249_0.1_results.csv"
main "$csv_file" "$num_rows" "$fastq_file"


# head -n 2001 AF237_0.1_poor_alignments.csv > AF237_0.1_poor_alignments.sampled.csv
# head -n 2001 AF237_0.1_results.csv > AF237_0.1_results.sampled.csv
# head -n 2001 AF249_0.1_poor_alignments.csv > AF249_0.1_poor_alignments.sampled.csv
# head -n 2001 AF249_0.1_results.csv > AF249_0.1_results.sampled.csv
