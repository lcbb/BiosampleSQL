import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve,auc
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import os

plt.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def roc_plot_from_counts(seq_counts_file, truth_table_file):
    # Load the seq_counts data from a CSV file
    seq_counts_df = pd.read_csv(seq_counts_file)
    # Convert the 'Match Count' column to a numpy array
    seq_counts = seq_counts_df['Match Count'].values

    # Load the truth_set data from a CSV file
    truth_set_df = pd.read_csv(truth_table_file)
    # Convert the second column to a numpy array
    truth_set = truth_set_df['Truth Value'].values

    # Create a binary classification based on your threshold
    #binary_classification = np.where(seq_counts >= threshold, 1, 0)

    # Compute ROC curve and ROC area
    fpr, tpr, _ = roc_curve(truth_set, seq_counts)
    roc_auc = auc(fpr, tpr)

    # Plot ROC curve
    plt.figure()
    lw = 2  # Line width
    plt.plot(fpr, tpr, color='darkorange', lw=lw, label='AUC ROC = %0.2f' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend(loc="lower right")

    # Extract the filename without extension to use as prefix for the figure
    output_prefix = seq_counts_file.split('.')[0]

    # Save the figure
    plt.savefig(f'{output_prefix}_roc_curve.png')

def roc_plot_from_multiple_counts(seq_counts_files, truth_table_files, combined_plot=True, figsize=(4, 4)):
    assert len(seq_counts_files) == len(truth_table_files), "Number of sequence count files should match with number of truth table files."

    lw = 1  # Line width

    if combined_plot:
        plt.figure(figsize=figsize)

    for seq_counts_file, truth_table_file in zip(seq_counts_files, truth_table_files):
        seq_counts_df = pd.read_csv(seq_counts_file)
        seq_counts = seq_counts_df['Match Count'].values

        truth_set_df = pd.read_csv(truth_table_file)
        truth_set = truth_set_df['Truth Value'].values

        fpr, tpr, _ = roc_curve(truth_set, seq_counts)
        roc_auc = auc(fpr, tpr)

        if not combined_plot:
            plt.figure(figsize=figsize)
        plt.plot([0, 1], [0, 1], color=[0.7, 0.7, 0.7], lw=lw, linestyle=':')
        plt.plot(fpr, tpr, color='darkorange', lw=lw, label=f'{os.path.splitext(os.path.basename(seq_counts_file))[0]} (AUC = {roc_auc:.2f})')
        

        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xticks([0, 0.5, 1.0])
        plt.yticks([0, 0.5, 1.0])
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.legend(loc="lower right")

        if not combined_plot:
            plt.savefig(f'{os.path.splitext(os.path.basename(seq_counts_file))[0]}_roc_curve.png', dpi=300)
            plt.close()

    if combined_plot:
        plt.savefig('combined_roc_curve.pdf')
        plt.close()

# Call the function
seq_counts_files = [
    r'match_analysis/2023_08_03_symptomatic_1_S2.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_age1_1_S3.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_age2_1_S4.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_age3_1_S5.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_symptomatic_2_S6.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_age1_2_S7.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_age2_2_S8.assembled_match_counts.csv',
    r'match_analysis/2023_08_03_age3_2_S9.assembled_match_counts.csv'
]

truth_table_files = [
    r'symp_neg_truthtable.csv',
    r'age1_truthtable.csv',
    r'age2_truthtable.csv',
    r'age3_truthtable.csv',
    r'symp_neg_truthtable.csv',
    r'age1_truthtable.csv',
    r'age2_truthtable.csv',
    r'age3_truthtable.csv'
]

roc_plot_from_multiple_counts(seq_counts_files, truth_table_files, combined_plot=False, figsize=(4, 4))