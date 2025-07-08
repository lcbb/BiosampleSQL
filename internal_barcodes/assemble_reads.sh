#!/bin/bash

# Set the input and output directories
INPUT_DIR=raw_reads
OUTPUT_DIR=assembled_reads
DEDUPE_DIR=deduplicated_reads

mkdir -p $OUTPUT_DIR
mkdir -p $DEDUPE_DIR

# Set the PEAR options
# MIN_OVERLAP=20
THREADS=20

# Loop over the input files
for file1 in $INPUT_DIR/*_R1_001.fastq.gz
do
  # Set the input file names
  file2=${file1/_R1/_R2}
  sample_name=$(basename ${file1} _R1_001.fastq.gz)

  # Set the output file name
  output_file=$OUTPUT_DIR/${sample_name}
  dedup_file=$DEDUPE_DIR/${sample_name}

  # Decompress the FastQ files
  gunzip -c $file1 > ${file1%.gz}
  gunzip -c $file2 > ${file2%.gz}

  # Let's count the total reads from the decompressed file
  total_reads=$(($(cat ${file1%.gz} | wc -l)/4))

  echo -e "${file1%.gz}\n${file2%.gz}" > input.fofn
  FastUniq -i input.fofn -t q -o ${dedup_file}_R1.fastq -p ${dedup_file}_R2.fastq
  rm input.fofn ${file1%.gz} ${file2%.gz}

  # Let's count the deduplicated reads
  deduplicated_reads=$(($(cat ${dedup_file}_R1.fastq | wc -l)/4))

  echo "Sample $sample_name: Total reads = $total_reads, Deduplicated reads = $deduplicated_reads"

  # Now, assemble the paired-end reads using PEAR
  # pear -f ${dedup_file}_R1.fastq -r ${dedup_file}_R2.fastq -o $output_file -v $MIN_OVERLAP -j $THREADS
  pear -f ${dedup_file}_R1.fastq -r ${dedup_file}_R2.fastq -o $output_file -j $THREADS

  # Optional: remove the unmerged files to save space
  rm $output_file.unassembled.forward.fastq $output_file.unassembled.reverse.fastq $output_file.discarded.fastq
done
