import argparse
import pandas as pd
import os
import sys

# Reverse complement function without using Bio library
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def extract_sequences(fastq_file, adapter_left, adapter_right):
    # Run cutadapt to extract only the region between adapters
    cutadapt_output = fastq_file.replace(".merged.fastq", ".cutadapt.fastq")
    command = f"""
    cutadapt -g '{adapter_left};min_overlap=9...{adapter_right};min_overlap=9' --rename '{{id}}_{{match_sequence}}_{{rc}}' --rc -e 0.1 --cores 8 --discard-untrimmed -o {cutadapt_output} {fastq_file}
    """
    
    # Run the command using os.system and check if it executes successfully
    if os.system(command) != 0:
        sys.exit("Error: cutadapt command failed. Please check that cutadapt is installed and your parameters are correct.")

    # Parse cutadapt output and process
    sequences, failed_reads = [], []
    with open(cutadapt_output, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq_str = f.readline().strip()  # The extracted sequence between adapters
            f.readline()  # Plus line
            f.readline()  # Quality line

            # Split into barcode (first 15 bp) and target (remaining sequence)
            barcode, target = seq_str[:15], seq_str[15:]

            # Filter based on barcode pattern and target length
            if (
                len(target) in range(21, 26) and
                all((barcode[i] in 'AT' if c == 'W' else True) for i, c in enumerate("NNWNNWNNWNNWNNW"))
            ):
                sequences.append({'barcode': barcode, 'target': target})
            else:
                failed_reads.append({'barcode': barcode, 'target': target})  # Save the failed ones.

    return sequences, failed_reads

def save_results(fastq_file, sequences, failed_reads):
    # Create a DataFrame and save as CSV
    results_df = pd.DataFrame(sequences)
    results_csv = fastq_file.replace(".merged.fastq", "__results.csv")
    results_df.to_csv(results_csv, index=False)

    # Save failed read IDs as a separate CSV
    failed_df = pd.DataFrame({'read_id': failed_reads})
    failed_csv = fastq_file.replace(".merged.fastq", "__failed_reads.csv")
    failed_df.to_csv(failed_csv, index=False)

def main():
    parser = argparse.ArgumentParser(description="Extract and filter sequences using adapters.")
    parser.add_argument('-f', '--fastq_file', required=True, help="Path to the input fastq file.")
    parser.add_argument('-aL', '--adapter_left', required=True, help="Left adapter sequence (10 bp).")
    parser.add_argument('-aR', '--adapter_right', required=True, help="Right adapter sequence (10 bp).")
    
    args = parser.parse_args()

    sequences, failed_reads = extract_sequences(args.fastq_file, args.adapter_left, args.adapter_right)
    save_results(args.fastq_file, sequences, failed_reads)

if __name__ == "__main__":
    main()
