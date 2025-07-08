#!/bin/bash

if [[ ! -f *.fai ]]; then
        bwa index NC_045512.2.fa
        samtools faidx NC_045512.2.fa
        samtools dict NC_045512.2.fa -o NC_045512.2.dict
fi

fname="220078317"

mkdir -p "$fname"
cd "$fname"

if [[ ! -f *.sam ]]; then
        bwa mem -t 64 \
        -R '@RG\tID:"$fname"\tLB:"$fname"\tPL:ILLUMINA\tPM:NEXTSEQ\tSM:"$fname"' \
        ../NC_045512.2.fa \
        ../raw_reads/220078317_217_S1_L001_R1_001.fastq.gz  \
        ../raw_reads/220078317_217_S1_L001_R2_001.fastq.gz \
        > "$fname".aligned.sam
        bwa mem -t 64 \
        -R '@RG\tID:"${fname}.enc\tLB:"${fname}.enc\tPL:ILLUMINA\tPM:NEXTSEQ\tSM:"${fname}.enc' \
        ../NC_045512.2.fa \
        ../raw_reads/220078317_E_217_S6_L001_R1_001.fastq.gz \
        ../raw_reads/220078317_E_217_S6_L001_R2_001.fastq.gz \
        > "$fname".enc.aligned.sam
fi

if [[ ! -f *.sorted.dedup.bam ]]; then
        ~/bin/gatk-4.3.0.0/gatk MarkDuplicatesSpark \
        -I "$fname".aligned.sam \
        -M "$fname".dedup_metrics.txt \
        -O "$fname".sorted.dedup.bam \
        --remove-all-duplicates true --spark-master local[*]
        ~/bin/gatk-4.3.0.0/gatk MarkDuplicatesSpark \
        -I "$fname".enc.aligned.sam \
        -M "$fname".enc.dedup_metrics.txt \
        -O "$fname".enc.sorted.dedup.bam \
        --remove-all-duplicates true --spark-master local[*]
fi

if [[ ! -f *.raw_variants.vcf ]]; then
        ~/bin/gatk-4.3.0.0/gatk HaplotypeCaller --native-pair-hmm-threads 20 \
        -R ../NC_045512.2.fa \
        -I "$fname".sorted.dedup.bam \
        -O "$fname".raw_variants.vcf
        ~/bin/gatk-4.3.0.0/gatk HaplotypeCaller --native-pair-hmm-threads 20 \
        -R ../NC_045512.2.fa \
        -I "$fname".enc.sorted.dedup.bam \
        -O "$fname".enc.raw_variants.vcf

        ~/bin/gatk-4.3.0.0/gatk VariantFiltration \
        -V "$fname".raw_variants.vcf \
        -O "$fname".filtered_variants.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40"
        ~/bin/gatk-4.3.0.0/gatk VariantFiltration \
        -V "$fname".enc.raw_variants.vcf \
        -O "$fname".enc.filtered_variants.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40"

        ~/bin/gatk-4.3.0.0/gatk SelectVariants --exclude-filtered \
        -V "$fname".filtered_variants.vcf \
        -O "$fname".hard-filtered.vcf
        ~/bin/gatk-4.3.0.0/gatk SelectVariants --exclude-filtered \
        -V "$fname".enc.filtered_variants.vcf \
        -O "$fname".enc.hard-filtered.vcf

        ~/bin/gatk-4.3.0.0/gatk FastaAlternateReferenceMaker -R ../NC_045512.2.fa -V "$fname".hard-filtered.vcf \
        -O "$fname".consensus.fasta
        ~/bin/gatk-4.3.0.0/gatk FastaAlternateReferenceMaker -R ../NC_045512.2.fa -V "$fname".enc.hard-filtered.vcf \
        -O "$fname".enc.consensus.fasta
fi

if [[ ! -f *.alignment_metrics.txt ]]; then
        ~/bin/gatk-4.3.0.0/gatk CollectAlignmentSummaryMetrics -R ../NC_045512.2.fa -I "$fname".sorted.dedup.bam \
        -O "$fname".alignment_metrics.txt
        ~/bin/gatk-4.3.0.0/gatk CollectAlignmentSummaryMetrics -R ../NC_045512.2.fa -I "$fname".enc.sorted.dedup.bam \
        -O "$fname".enc.alignment_metrics.txt
fi

if [[ ! -f *.insert_metrics.txt ]]; then
        ~/bin/gatk-4.3.0.0/gatk CollectInsertSizeMetrics -I "$fname".sorted.dedup.bam \
        -O "$fname".insert_metrics.txt -H "$fname".insert_size_histogram.pdf
        ~/bin/gatk-4.3.0.0/gatk CollectInsertSizeMetrics -I "$fname".enc.sorted.dedup.bam \
        -O "$fname".enc.insert_metrics.txt -H "$fname".enc.insert_size_histogram.pdf
fi

if [[ ! -f *.coverage.txt ]]; then
        samtools ampliconstats ../neb_vss2a.primer.bed "$fname".sorted.dedup.bam > "$fname".coverage.txt
        samtools ampliconstats ../neb_vss2a.primer.bed "$fname".enc.sorted.dedup.bam > "$fname".enc.coverage.txt
fi

cd ..