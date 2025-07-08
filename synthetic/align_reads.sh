#!/bin/bash

FASTQ_DIR="01_reads"

mkdir -p 02_bam
mkdir -p 03_vcf
mkdir -p 04_freyja

for main_prefix in $(ls 01_reads/*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//g' | sort | uniq); do
    # Extract only the filename without path for saving in the right directory
    base_prefix=$(basename $main_prefix)
    
    # Create read group info
    rg="@RG\tID:${base_prefix}\tSM:${base_prefix}\tPL:Illumina"
    
    # Sequence alignment    
    bwa mem -t 64 -R $rg 00_ref/NC_045512.2.fa 01_reads/"${base_prefix}_R1_001.fastq.gz" 01_reads/"${base_prefix}_R2_001.fastq.gz" | samtools view -bS - | samtools sort -@ 64 -o 02_bam/"${base_prefix}_sorted.bam" -
    
    # GATK MarkDuplicates
    ~/bin/gatk-4.3.0.0/gatk MarkDuplicates \
    -I 02_bam/"${base_prefix}_sorted.bam" \
    -O 02_bam/"${base_prefix}_sorted_nodup.bam" \
    -M 02_bam/"${base_prefix}_dup_metrics.txt" \
    --CREATE_INDEX true
    
    # Use LoFreq for variant calling
    lofreq call-parallel --pp-threads 64 -f 00_ref/NC_045512.2.fa -o 03_vcf/"${base_prefix}.vcf" 02_bam/"${base_prefix}_sorted_nodup.bam"
    
    # Use samtools to calculate read depth
    samtools mpileup -aa -A -d 600000 -Q 10 -q 0 -B -f 00_ref/NC_045512.2.fa 02_bam/"${base_prefix}_sorted.bam" | cut -f1-4 > 04_freyja/"${base_prefix}_depths.txt"

    # Use freyja for variant calling and variant demix
    freyja demix 03_vcf/"${base_prefix}.vcf" 04_freyja/"${base_prefix}_depths.txt" --output 04_freyja/"${base_prefix}_demix_output.txt"
done

