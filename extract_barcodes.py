import gzip
import pandas as pd
from Bio import SeqIO
from skbio.alignment import StripedSmithWaterman
import argparse
from Bio.Seq import Seq
import re
import os

# Update! See main() for template sequences. 
# Define the template sequence
# template_sequence = "ACGCGAGCTGCATGTGTCAGANNWNNWNNWNNWNNWGGACAGTAAGAAGGAAAAACAGGCTTCTTCAAG"
# extended tamplate_sequence = "ACGCGAGCTGCATGTGTCAGANNWNNWNNWNNWNNWGGACTGAGGGCCATGGACACGGGCTTCTTCAAGTCCGCCATGCCCGA"

# Positions of the different parts within the template

universal_primer_start = 0
universal_primer_end = 21
barcode_start = 21
barcode_end = 36
target_start = 36
target_end = 59
flanking_start = 59
flanking_end = 83

def extract_parts(sequence, template, flanking_n=10, min_matches=9):
    # Reverse complement the sequence
    rev_comp_sequence = str(Seq(sequence).reverse_complement())

    # Align with original sequence
    aligner = StripedSmithWaterman(
        sequence,
        gap_open_penalty=10,
        gap_extend_penalty=2,
        match_score=2,
        mismatch_score=-1
    )
    alignment_orig = aligner(template)

    # Align with reverse complement sequence
    aligner_rev = StripedSmithWaterman(
        rev_comp_sequence,
        gap_open_penalty=10,
        gap_extend_penalty=2,
        match_score=2,
        mismatch_score=-1
    )
    alignment_rev = aligner_rev(template)

    # Choose the best alignment
    alignment = alignment_orig if alignment_orig.optimal_alignment_score > alignment_rev.optimal_alignment_score else alignment_rev

    # Aligned sequences from the best alignment
    aligned_sequence = alignment.aligned_query_sequence
    aligned_template = alignment.target_sequence
    cigar = alignment.cigar
    # print(alignment.target_begin)
    # print(aligned_sequence)
    # print(cigar)

    # Manually adjust the template and sequence to insert gaps based on the CIGAR string
    adjusted_template = []
    adjusted_sequence = []
    template_idx = 0
    sequence_idx = 0

    # Use regex to correctly parse the CIGAR string
    cigar_operations = re.findall(r'(\d+)([MID])', cigar)
    
    adjusted_sequence.extend(['-'] * alignment.target_begin)
    for length, op in cigar_operations:
        length = int(length)
        if op == 'M':  # Match or mismatch
            adjusted_template.extend(aligned_template[template_idx:template_idx+length])
            adjusted_sequence.extend(aligned_sequence[sequence_idx:sequence_idx+length])
            template_idx += length
            sequence_idx += length
        elif op == 'I':  # Insertion in the query
            adjusted_template.extend(['-'] * length)
            adjusted_sequence.extend(aligned_sequence[sequence_idx:sequence_idx+length])
            sequence_idx += length
        elif op == 'D':  # Deletion in the query
            adjusted_template.extend(aligned_template[template_idx:template_idx+length])
            adjusted_sequence.extend(aligned_sequence[sequence_idx:sequence_idx+length])
            # adjusted_sequence.extend(['-'] * length)
            sequence_idx += length
            template_idx += length

    # Convert lists back to strings
    adjusted_template = ''.join(adjusted_template)
    adjusted_sequence = ''.join(adjusted_sequence)
    # print("Adjusted Aligned Sequence: ", adjusted_sequence)
    # print("Adjusted Aligned Template: ", adjusted_template)
    
    # print(str(len(adjusted_template)))
    if (len(adjusted_template) < 58) or (len(adjusted_sequence) < 58):
        return None, alignment
    
    def extract_query_region(start, end):
        adjusted_start = start
        adjusted_end = end
        
        if start > len(adjusted_sequence):
            return ""
        elif (end > len(adjusted_sequence)) or (end > len(adjusted_template)):
            end = min([len(adjusted_sequence), len(adjusted_template)])
        
        # Adjust for gaps (insertions in the template or deletions in the query) before the start position
        for i in range(start):
            if adjusted_template[i] == '-':  # Gap in the template
                adjusted_start += 1  # Skip this position in the query
                adjusted_end += 1  # Extend the region to include the equivalent length in the query
            # elif aligned_sequence[i] == '-':  # Insertion in the query
                # adjusted_start += 1  # Skip this position in the template

        # Adjust for gaps within the region
        for i in range(start, end):
            if adjusted_template[i] == '-':  # Deletion in the template
                adjusted_end += 1  # Extend the region to include the equivalent length in the query
            # if aligned_sequence[i] == '-':  # Insertion in the query
            #     adjusted_end += 1  # Extend the region to include the equivalent length in the query

        # Extract and return the query sequence from the adjusted start to adjusted end
        return aligned_sequence[adjusted_start:adjusted_end]
    
    def extract_template_region(start, end):
        adjusted_start = start
        adjusted_end = end
        # print(f"length of adjusted template: {len(adjusted_template)}")
        # print(f"length of adjusted sequence: {len(adjusted_sequence)}")
        # print(f"start: {str(start)}")
        # print(f"end: {str(end)}")
        if start > len(adjusted_template):
            return ""
        if end > len(adjusted_sequence):
            end = len(adjusted_sequence) 
        if end > len(adjusted_template):
            end = len(adjusted_template)
        
        # Adjust for gaps (insertions in the query or deletions in the template) before the start position
        for i in range(start):
            # if aligned_sequence[i] == '-':  # Insertion in the query
                # adjusted_start += 1  # Skip this position in the template
                # adjusted_end += 1  # Extend the region to include the equivalent length in the template
            if adjusted_template[i] == '-':  # Deletion in the template
                adjusted_start += 1  # Skip this position in the template

        # Adjust for gaps within the region
        for i in range(start, end):
        #     if aligned_sequence[i] == '-':  # Insertion in the query
        #         adjusted_end += 1  # Extend the region to include the equivalent length in the template
            if adjusted_template[i] == '-':  # Deletion in the template
                adjusted_end += 1  # Extend the region to include the equivalent length in the template

        # Extract and return the template sequence from the adjusted start to adjusted end
        return adjusted_template[adjusted_start:adjusted_end]


    x = alignment.target_begin
    
    universal_primer = extract_query_region(universal_primer_start, universal_primer_end-x)
    barcode = extract_query_region(barcode_start-x, barcode_end-x)
    target = extract_query_region(target_start-x, target_end-x)
    flanking = extract_query_region(flanking_start-x, flanking_end-x)
    
    template_flanking = extract_template_region(flanking_start, flanking_end)
    template_primer = extract_template_region(universal_primer_start, universal_primer_end)
    
    # print(universal_primer[-10:], template_primer[-10:])
    # print(flanking[0:10], template_flanking[0:10])
    
    if (len(universal_primer) < 10) or (len(flanking) < 10):
        return None, alignment
    
    matches_primer = sum(1 for a, b in zip(universal_primer[-10:], template_primer[-10:]) if a == b)
    matches_flanking = sum(1 for a, b in zip(flanking[0:10], template_flanking[0:10]) if a == b)
    
    
    if matches_primer < min_matches or matches_flanking < min_matches:
        return None, alignment
    
    
    return {
        'universal_primer': universal_primer,
        'barcode': barcode,
        'target': target,
        'flanking': flanking
    }, alignment



# Read and process each FASTQ file
def process_fastq(file_path, template):
    results = []
    poor_alignments = []
    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequence = str(record.seq)
            parts, alignment = extract_parts(sequence, template)
            if parts:
                parts.update({
                    'read_id': record.id,
                    'average_quality': sum(record.letter_annotations["phred_quality"]) / len(record),
                    'alignment_score': alignment.optimal_alignment_score,
                    'read_length': len(sequence),
                    'sequence': sequence,
                    'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
                })
                results.append(parts)
            else:
                poor_alignments.append({
                    'read_id': record.id,
                    'average_quality': sum(record.letter_annotations["phred_quality"]) / len(record),
                    'alignment_score': alignment.optimal_alignment_score,
                    'read_length': len(sequence),
                    'sequence': sequence
                })
    return results, poor_alignments

def save_to_csv(results, poor_alignments, output_prefix):
    df_results = pd.DataFrame(results)
    df_poor_alignments = pd.DataFrame(poor_alignments)
    df_results.to_csv(f"results/{output_prefix}_results.csv", index=False)
    df_poor_alignments.to_csv(f"results/{output_prefix}_poor_alignments.csv", index=False)

def main():
    parser = argparse.ArgumentParser(description='Process FASTQ files for Smith-Waterman alignment.')
    parser.add_argument('fastq_files', nargs='+', help='Paths to the FASTQ files to process.')
    args = parser.parse_args()
    
    templates = pd.read_csv("template_sequences.csv")
    os.makedirs("results", exist_ok=True)
    for file in args.fastq_files:
        print(f"Processing {file}...")
        sample_id = file.split('/')[-1].split('.')[0]
        filtered_templates = templates[templates["Sample"] == sample_id].reset_index()
        
        if not filtered_templates.empty:
            template_sequence = filtered_templates["Template"].iloc[0]
            # Proceed with the template_sequence
            print(f"Template sequence for {sample_id}: {template_sequence}")
        else:
            print(f"No template sequence found for sample {sample_id}")
        
        results, poor_alignments = process_fastq(file, template_sequence)
        print(f"Saving {file}...")
        output_name = file.split('/')[-1]  # Extracts the filename from the full path
        output_name = "_".join([output_name.split('.')[0], ".".join(output_name.split('.')[3:5])])

        save_to_csv(results, poor_alignments, output_name)
        print("Done.")

if __name__ == "__main__":
    main()
    

# module load conda3/202311
# python extract_barcodes.py ../070724_nextseq_DP_merged/8e_rep1_S5.extendedFrags.fastq

# bsub -q standard -e run_stats/extract_seqtik.error -o run_stats/extract_seqtik.out -M 200000MB -n 1 python extract_barcodes.py merged_fastq_seqtik/*.merged.fastq 

