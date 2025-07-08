#!/bin/bash

mkdir -p 05_stats

while read -r bamfile; do
    filename=$(basename -- "$bamfile")
    filename_no_ext="${filename%.*}"
    
    samtools ampliconstats nCoV-2019.primer.bed 02_bam/$bamfile > 05_stats/"${filename_no_ext}_ampliconstats.txt"
    samtools depth -a 02_bam/$bamfile > 05_stats/"${filename_no_ext}_depth.txt"
    
    echo "Processed $bamfile"
    
 done < "bamfiles.txt"
