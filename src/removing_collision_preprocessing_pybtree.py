import os
import glob
import pandas as pd
import numpy as np
import pandas as pd
import ast  # Use ast.literal_eval to safely evaluate string representations of lists
from pybktree import BKTree
from polyleven import levenshtein

# Run an exploratory_analysis_testing file first to get the results (not combined). Using this results, we will filter them

# Function to update the DataFrame based on Levenshtein distance filtering
def filter_Cas9_unique_targets(df):
    for index, row in df.iterrows():
        control_seq = row['control_target_seq']  # Control sequence
        unique_targets = row['Cas9_unique_targets']  # List of Cas9 targets
        counts_per_target = row['Cas9_counts_per_target']  # Corresponding counts
        
        # Ensure Cas9_unique_targets and Cas9_counts_per_target are treated as lists
        if isinstance(unique_targets, str):
            unique_targets = ast.literal_eval(unique_targets)
        if isinstance(counts_per_target, str):
            counts_per_target = ast.literal_eval(counts_per_target)
        
        # Handle case where counts_per_target might be a float or other non-iterable
        if isinstance(counts_per_target, (float, int)):
            counts_per_target = [counts_per_target]  # Convert to list if it's a single value
        
        # Handle empty or invalid unique_targets (skip row if necessary)
        if not isinstance(unique_targets, list) or len(unique_targets) == 0:
            # Skip processing this row if unique_targets is not a valid list
            df.at[index, 'Cas9_unique_targets'] = []
            df.at[index, 'Cas9_counts_per_target'] = []
            df.at[index, 'Cas9_n_targets'] = 0
            df.at[index, 'Cas9_counts'] = 0
            continue  # Skip to the next row
        
        # Lists to store filtered targets and counts
        filtered_targets = []
        filtered_counts = []
        
        # Calculate edit distances and filter based on a distance threshold of 2
        for i, target in enumerate(unique_targets):
            # Calculate the edit distance between control_seq and each Cas9 target
            edit_dist = levenshtein(control_seq, target)
            
            if edit_dist <= 2:  # Keep only if edit distance is <= 2
                filtered_targets.append(target)
                filtered_counts.append(counts_per_target[i])
        
        # Update the DataFrame with the filtered results
        df.at[index, 'Cas9_unique_targets'] = filtered_targets
        df.at[index, 'Cas9_counts_per_target'] = filtered_counts
        df.at[index, 'Cas9_n_targets'] = len(filtered_targets)
        df.at[index, 'Cas9_counts'] = sum(filtered_counts)  # Ensure counts are summed as integers

    return df

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



def list_files_not_matching_pattern(folder_path, pattern):
    # Get all files in the folder
    all_files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    
    # Get files matching the pattern
    matching_files = set(os.path.basename(f) for f in glob.glob(os.path.join(folder_path, pattern)))
    
    # Filter out matching files
    non_matching_files = [f for f in all_files if f not in matching_files]
    
    return non_matching_files

# Usage
folder_path = 'results_log2FC'  # Replace with your folder path
pattern = '*_combined.csv'  # Replace with the specific pattern, e.g., 'file_*_2024.txt'
sample_results_files = list_files_not_matching_pattern(folder_path, pattern)


log2FC_pseudocount = 10

for result_file in sample_results_files:
    output_filename = result_file.replace('.csv', '')
    print(f'Processing {output_filename}.')
    result_df = pd.read_csv(f'{folder_path}/{result_file}')
    filtered_result_df = filter_Cas9_unique_targets(result_df)
    filtered_results_combined_df = combine_same_targets(filtered_result_df)
    
    filtered_result_df.to_csv(f"{folder_path}/{output_filename}_filtered.csv", index=False)
    filtered_results_combined_df.to_csv(f"{folder_path}/{output_filename}_filtered_combined.csv", index=False)
    

