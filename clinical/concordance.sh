#!/bin/bash

samples="220078317 22MG14803 22MG14761 22MG14613 22MG4776"

# Initialize the output file with headers
echo -e "Sample\tPrecision\tRecall\tF1_Score" > concordance_results.tsv

for sample in $samples; do
    echo "Processing sample $sample..."  # Check 1: Is the script entering the loop?
    
    cd "$sample"
    
    if [[ -r "$sample.enc.hard-filtered.vcf" && -r "$sample.hard-filtered.vcf" ]]; then  # Check 2: Do the input files exist and are they readable?
        # Normalize and decompose both VCFs using BCFtools
        bcftools norm -m -both $sample.enc.hard-filtered.vcf -Oz -o $sample.enc.hard-filtered.norm.vcf.gz
        echo "bcftools norm exit status: $?"  # Check 3: Did bcftools norm succeed?
        
        bcftools norm -m -both $sample.hard-filtered.vcf -Oz -o $sample.hard-filtered.norm.vcf.gz
        echo "bcftools norm exit status: $?"  # Check 3: Did bcftools norm succeed?
        
        # Generate index files
        bcftools index $sample.enc.hard-filtered.norm.vcf.gz
        echo "bcftools index exit status: $?"  # Check 3: Did bcftools index succeed?
        
        bcftools index $sample.hard-filtered.norm.vcf.gz
        echo "bcftools index exit status: $?"  # Check 3: Did bcftools index succeed?
        
        # Compare the VCFs using BCFtools
        bcftools isec -n=2 -w1 -Oz -o $sample.overlap.vcf.gz $sample.enc.hard-filtered.norm.vcf.gz $sample.hard-filtered.norm.vcf.gz
        echo "bcftools isec exit status: $?"  # Check 3: Did bcftools isec succeed?
        
        # Count the number of variants in each file
        TP=$(bcftools view -H $sample.overlap.vcf.gz | wc -l)
        FP=$(bcftools view -H $sample.enc.hard-filtered.norm.vcf.gz | wc -l)
        FN=$(bcftools view -H $sample.hard-filtered.norm.vcf.gz | wc -l)
        
        # Subtract the TP from the FP and FN to get the correct numbers
        FP=$((FP-TP))
        FN=$((FN-TP))
        
        echo "TP=$TP, FP=$FP, FN=$FN"  # Check 4: What are the values of TP, FP, and FN?
        
        # Now calculate precision, recall, and F1 score
        precision=$(echo "$TP / ($TP + $FP)" | bc -l)
        recall=$(echo "$TP / ($TP + $FN)" | bc -l)
        f1_score=$(echo "2 * (($precision * $recall) / ($precision + $recall))" | bc -l)
        
        # Write the results to the output file
        echo -e "$sample\t$precision\t$recall\t$f1_score" >> ../concordance_results.tsv
    else
        echo "Input files for $sample not found or not readable, skipping..."
    fi

    cd ..
done
