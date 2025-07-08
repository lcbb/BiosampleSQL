#!/bin/bash

# Define the directory to search
DIRECTORY="05_stats"

# Define the file to save the output to
OUTPUT_FILE="samtools.grep.tsv"

# Grep all files in the directory and save the output to a file
grep -r "^FPCOV" "$DIRECTORY" > "$OUTPUT_FILE"

# Transpose the file using awk
awk '{ for(i=1;i<=NF;i++) a[i,NR]=$i } NF>p { p=NF } END { for(j=1;j<=p;j++) { str=a[j,1]; for(i=2;i<=NR;i++) str=str"\t"a[j,i]; print str } }' "$OUTPUT_FILE" > "${DIRECTORY}.${OUTPUT_FILE}.transpose.tsv"

# Clean up the original output file
rm "$OUTPUT_FILE"
