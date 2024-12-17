import argparse
import pandas as pd
import numpy as np
import os
import math

# Generate a results.csv table

df_combined = pd.read_csv(f'results_cutadaptv2/all_results_combined_cleaned_21_25.csv')

def split_df(df, label):
    # Filter for the specific label and treatment
    df_filtered = df[(df['label'] == label)]
    return df_filtered

labels = df_combined['label'].unique().tolist()
print(f'labels: {labels}')
indiv_dfs = {}
for label in labels:
    indiv_dfs[label] = split_df(df_combined, label)

# Define minimum read count
minimum_read_count = 1
log2FC_pseudocount = 1
save_folder="results_log2FC"

# Function to get barcode-target counts for good barcodes
def get_good_barcodes(df):
    
    # Group by barcode and target, and count occurrences
    barcode_target_counts = df.groupby(['barcode', 'target']).size().reset_index(name='count')
    
    # Find barcodes that are associated with only one unique target
    barcode_unique_target_counts = barcode_target_counts.groupby('barcode')['target'].nunique().reset_index(name='unique_target_count')
    good_barcodes = barcode_unique_target_counts[barcode_unique_target_counts['unique_target_count'] == 1]['barcode']
    
    # Filter the original counts to only include good barcodes
    barcode_target_filtered = barcode_target_counts[barcode_target_counts['barcode'].isin(good_barcodes)]
    
    # Further filter to keep barcode-target combinations with count >= minimum_read_count
    barcode_target_filtered = barcode_target_filtered[barcode_target_filtered['count'] >= minimum_read_count]
    
    # Return the filtered barcode-target counts
    return barcode_target_filtered


def count_mismatches(seq, target):
    mismatch_count = 0
    for s_base, t_base in zip(seq, target):
        if t_base == 'N':
            continue  # N can match any base
        elif s_base != t_base:
            mismatch_count += 1
    return mismatch_count


def cas9_with_good_barcodes(df, good_barcodes_df, new_target_seq):
    
    good_barcodes = good_barcodes_df['barcode'].unique()
    
    # Only focus on the good barcodes from control
    df_filtered = df[df['barcode'].isin(good_barcodes)]
    
    # Group by barcode and target to count occurrences of each target per barcode
    target_counts = df_filtered.groupby(['barcode', 'target']).size().reset_index(name='target_count')
    
    # Group by barcode to aggregate unique targets, number of unique targets, total row count, and per-target row count
    barcode_target_counts = target_counts.groupby('barcode').agg(
        Cas9_unique_targets=('target', lambda x: list(x)),  # List of unique targets
        Cas9_n_targets=('target', 'nunique'),  # Count of unique targets
        Cas9_counts_per_target=('target_count', lambda x: list(x)),  # List of row counts per target
        Cas9_counts=('target_count', 'sum')  # Total number of rows (occurrences) for each barcode
    ).reset_index()
    
    # Merge with good_barcodes_df
    results = pd.merge(good_barcodes_df, barcode_target_counts, on='barcode', how='left')
    
    # Rename the target and count columns
    results = results.rename(columns={'target': 'control_target_seq', 'count': 'control_counts'})
    
    results['control_ratios'] = results['control_counts']/results['control_counts'].sum()
    results['Cas9_ratios'] = results['Cas9_counts']/results['Cas9_counts'].sum()
    
    results['control_norm'] = results['control_ratios']*1000000
    results['Cas9_norm'] = results['Cas9_ratios']*1000000
    
    # results['log2FC'] = np.log2((results['Cas9_norm']+log2FC_pseudocount) / (results['control_norm']+log2FC_pseudocount)) ############################
    # results = results.sort_values(by='log2FC', ascending=False)
    
    results = results[results['control_target_seq'].str.len() == 23] 
    
    results['MM'] = results['control_target_seq'].apply(lambda x: count_mismatches(x, new_target_seq))
    
    return results

def combine_same_targets(df_results):
    # Group by control_target_seq and sum the control_counts and Cas9_counts
    df_results_combined = df_results.groupby('control_target_seq').agg(
        control_counts=('control_counts', 'sum'),
        Cas9_counts=('Cas9_counts', 'sum'),
        MM=('MM', 'first')
    ).reset_index()
    
    df_results_combined['control_ratios'] = df_results_combined['control_counts']/df_results_combined['control_counts'].sum()
    df_results_combined['Cas9_ratios'] = df_results_combined['Cas9_counts']/df_results_combined['Cas9_counts'].sum()
    
    df_results_combined['control_norm'] = df_results_combined['control_ratios']*1000000
    df_results_combined['Cas9_norm'] = df_results_combined['Cas9_ratios']*1000000
    
    # Calculate log2FC for log2(normalized_Cas9_counts / normalized_control_counts)
    df_results_combined['log2FC'] = np.log2((df_results_combined['Cas9_norm']+ log2FC_pseudocount) / (df_results_combined['control_norm']+log2FC_pseudocount))
    df_results_combined = df_results_combined.sort_values('log2FC', ascending=False)
    return df_results_combined

def off_target_analysis_pipeline(df_control, df_treated, name):
    target_seq = df_control['Sequence'].iloc[0]
    new_target_seq = target_seq[:20] + 'N' + target_seq[21:]
    
    good_barcodes_control = get_good_barcodes(df_control)
    
    df_results = cas9_with_good_barcodes(df_treated, good_barcodes_control, new_target_seq)
    df_results_combined = combine_same_targets(df_results)
    
    df_results.to_csv(f"{save_folder}/results_{name}.csv", index=False)
    df_results_combined.to_csv(f"{save_folder}/results_{name}_combined.csv", index=False)
    
    # return df_results, df_results_combined



unique_labels = combined_df['Label'].unique()

for label in unique_labels:
    df_control = combined_df[(combined_df['Label'] == label) & (combined_df['Treatment'] == 'SrfI')]
    df_treated = combined_df[(combined_df['Label'] == label) & (combined_df['Treatment'] == 'Cas9')]
    
    # Only run if both dataframes are non-empty
    if not df_control.empty and not df_treated.empty:
        off_target_analysis_pipeline(df_control, df_treated, name=label)