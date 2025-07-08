import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def create_heatmaps(match_file_path, truth_file_path, scale_type='linear'):
    # Read the CSV file
    df_match = pd.read_csv(match_file_path)
    df_truth = pd.read_csv(truth_file_path)

    # Remove 'barcode' from the Sequence ID
    df_match['Sequence ID'] = df_match['Sequence ID'].str.replace('barcode', '').astype(int)
    df_truth['Sequence ID'] = df_truth['Sequence ID'].str.replace('barcode', '').astype(int)

    # Generate a complete range of Sequence IDs and merge with the existing DataFrame
    all_ids = pd.DataFrame({'Sequence ID': range(0, 96)})
    df_match = all_ids.merge(df_match, on='Sequence ID', how='left').fillna(0)
    df_truth = all_ids.merge(df_truth, on='Sequence ID', how='left').fillna(0)

    # Create Sequence ID grid
    seq_id_grid = np.arange(96).reshape(8, 12)

    # Prepare the data for the heatmaps
    match_heatmap_data = pd.pivot_table(df_match, values='Match Count', index=(df_match['Sequence ID'] // 12), columns=(df_match['Sequence ID'] % 12), fill_value=0)
    truth_heatmap_data = pd.pivot_table(df_truth, values='Truth Value', index=(df_truth['Sequence ID'] // 12), columns=(df_truth['Sequence ID'] % 12), fill_value=0)

    # Apply the scaling
    if scale_type == 'log10':
        match_heatmap_data = np.log10(match_heatmap_data + 1)  # +1 to avoid division by zero
    elif scale_type == 'log2':
        match_heatmap_data = np.log2(match_heatmap_data + 1)  # +1 to avoid division by zero
    # For linear scale, do nothing

    # Create the heatmaps
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4))

    sns.heatmap(truth_heatmap_data, cmap='viridis', ax=ax1, cbar_kws={'label': 'Truth value'}, annot=seq_id_grid, fmt="d", linewidths=1.0, linecolor="white")
    ax1.set_title('Expected')

    sns.heatmap(match_heatmap_data, cmap='viridis', ax=ax2, cbar_kws={'label': 'Match counts ({})'.format(scale_type)}, annot=seq_id_grid, fmt="d", linewidths=1.0, linecolor="white")
    ax2.set_title('Experiment')

    for ax in (ax1, ax2):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')

    # Save the heatmap as a PNG file
    file_name = os.path.splitext(os.path.basename(match_file_path))[0]
    output_file = f'{file_name}_heatmap_{scale_type}.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()


# Specify the filenames
match_files = [
    r'match_analysis/2023_08_02_symptomatic_S2.assembled_match_counts.csv',
    r'match_analysis/2023_08_02_age1_S3.assembled_match_counts.csv',
    r'match_analysis/2023_08_02_age2_S4.assembled_match_counts.csv',
    r'match_analysis/2023_08_02_age3_S5.assembled_match_counts.csv',
]

truth_files = [
    r'symp_neg_truthtable.csv',
    r'age1_truthtable.csv',
    r'age2_truthtable.csv',
    r'age3_truthtable.csv'
]

for match_file, truth_file in zip(match_files, truth_files):
    create_heatmaps(match_file, truth_file, scale_type='log2')