from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import csv
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import re

# Function to sort alphanumerically
def natural_sort(s):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]

# Input and output directories
INPUT_DIR = "assembled_reads"
OUTPUT_DIR = "match_analysis"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Number of reads to parse
max_reads = 1_000_000

# Load the FASTA file into a dictionary
fa = "templates_full.fasta"
fasta_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

# Create a PairwiseAligner object
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1
aligner.mismatch_score = -20
aligner.open_gap_score = -2
aligner.extend_gap_score = -2

# Loop over the files in the input directory
for filename in os.listdir(INPUT_DIR):

    if not filename.endswith(".fastq") and not filename.endswith(".fq"):
        continue

    # Set the input and output filenames
    input_file = os.path.join(INPUT_DIR, filename)
    output_file = os.path.join(OUTPUT_DIR, os.path.splitext(filename)[0])

    # Loop over the reads and count the number of matches for each sequence in the FASTA file
    match_counts = defaultdict(int, {k: 0 for k in fasta_dict.keys()})

    with open(input_file, "r") as handle:
        pbar = tqdm(total=max_reads, desc="Processing " + input_file)
        for record in SeqIO.parse(handle, "fastq"):

            # Calculate the alignment scores for this read against each template sequence
            alignment_scores = [aligner.score(record.seq, fasta_seq) for fasta_seq in fasta_dict.values()]

            # Set the alignment score threshold to consider a match
            threshold = 35

            # If the alignment score is above the threshold, count it as a match
            for i, fasta_id in enumerate(fasta_dict.keys()):
                if alignment_scores[i] >= threshold:
                    match_counts[fasta_id] += 1

            if max_reads is not None:
                pbar.update(1)
                if pbar.n >= max_reads:
                    break
        pbar.close()

    # Output the match counts to a CSV file
    with open(os.path.join(OUTPUT_DIR, os.path.splitext(filename)[0] + "_match_counts_1M.csv"), "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Sequence ID", "Match Count"])
        for seq_id, count in sorted(match_counts.items(), key=lambda x: natural_sort(x[0])):
            writer.writerow([seq_id, count])

    # Plot counts
    fid = []
    cnt = []

    for fasta_id in fasta_dict:
       fid.append(fasta_id)
       cnt.append(match_counts[fasta_id])

    # Plot match counts
    plt.figure(figsize=(24,2))

    ylims = (0, max(cnt)*1.5)

    plt.bar(range(len(match_counts)), cnt)
    plt.xticks(range(len(match_counts)), fid, rotation=90, fontsize=10)

    for i, c in enumerate(cnt):
        plt.text(i, ylims[1]/50, str(c), horizontalalignment="center", verticalalignment="bottom", rotation="vertical")

    plt.ylim(ylims)

    plt.savefig(output_file + "_counts_1M" + ".png", dpi=300, bbox_inches="tight")
    plt.close()