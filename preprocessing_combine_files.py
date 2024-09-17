import pandas as pd
import os
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process folder path input.")
parser.add_argument('--folder_path', type=str, required=True, help='Folder path where results are stored')

# Parse the arguments
args = parser.parse_args()

# Load sample information
sample_info = pd.read_csv('sample_info.csv')

# Use the input folder path
folder_path = args.folder_path

run = True

if run == True:
    # List to hold individual DataFrames
    data_frames = []
    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith("_results.csv"):
            # Read the CSV file
            print(f"Reading: {filename}")
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path)

            # Extract the sample name from the filename (assuming format: AF211_S1_results.csv)
            sample_name = filename.split('_')[0]  # Change this if the sample name extraction rule is different
            sampling_ratio = filename.split('_')[1]

            # Add the "Sample" column
            df['sample'] = sample_name
            df['sampling_ratio'] = sampling_ratio
            df = df[['barcode','target','sample', 'sampling_ratio']]
            # Append the DataFrame to the list
            data_frames.append(df)
            

    # Concatenate all DataFrames
    df_combined = pd.concat(data_frames, ignore_index=True)
    
    # sample_info = pd.read_csv("template_sequences.csv")
    sample_info = pd.read_csv("sample_info.csv")

    # Merge the combined DataFrame with the lookup table
    df_combined = pd.merge(df_combined, sample_info, left_on='sample', right_on = "Sample", how='left')

    df_combined.drop(columns=['Sample'])
    df_combined = df_combined[['barcode', 'target', 'Sequence', 'sample', 'sampling_ratio', 'Treatment', 'Group', 'Target', 'Concentration', 'UMI_format']]
    
    print("Saving combined file...")
    df_combined.to_csv(f'{folder_path}/all_results_combined.csv', index=False) # Export as a slightly less big file for downstream analysis.

print("Cleaning up and filtering combined file...")
df_combined_clean = df_combined.copy()

# Clean up some NAs and replace - with ''.
df_combined_clean = df_combined_clean.dropna()
df_combined_clean = df_combined_clean.replace('-', '', regex=True)

# Filter for only barcodes that are 15 bp and randomized targets that are between 21 and 30 bp
# Remove reads with sequence ambiguity (Ns in the reads)
df_combined_clean['barcode'] = df_combined_clean['barcode'].astype(str)
df_combined_clean['target'] = df_combined_clean['target'].astype(str)
df_combined_clean = df_combined_clean[~df_combined_clean['barcode'].str.contains('N') & ~df_combined_clean['target'].str.contains('N')]
df_combined_clean = df_combined_clean[(df_combined_clean['barcode'].str.len()==15) & (df_combined_clean['target'].str.len()>=21) & (df_combined_clean['target'].str.len()<=30)]

# Add a label for important features
df_combined_clean['label'] = df_combined_clean['Target'] + '_' + df_combined_clean['Concentration'].astype(str) + '_' + df_combined_clean['UMI_format']

# Save results
print("Saving filterd combined file as all_results_combined_new_15_21-30_noNs.csv...")
df_combined_clean.to_csv(f'{folder_path}/all_results_combined_new_15_21-30_noNs.csv', index=False)
print("Done.")



# Unique read
print("Final touches: creating summary table for plotting...")


df_combined_clean = pd.read_csv(f'{folder_path}/all_results_combined_new_15_21-30_noNs.csv')

df_combined_clean['unique_read'] = df_combined_clean['barcode']+df_combined_clean['target']

# Group by label and sampling_ratio, then count the unique reads and total reads
summary_df = df_combined_clean.groupby(['label', 'sampling_ratio']).agg(
    unique_read_count=('unique_read', 'nunique'),
    total_read_count=('unique_read', 'count')
).reset_index()

# Display the result
summary_df.head()
summary_df.to_csv(f'{folder_path}/summary_counts_for_plot.csv', index=False)

print('Done: summary_counts_for_plot.csv')