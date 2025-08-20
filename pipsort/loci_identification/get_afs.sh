#!/bin/bash

# File containing list of vid strings
VIDS_FILE="vids.txt"
INPUT_FILE="vat_complete_v7.1.bgz.tsv.gz"
OUTPUT_FILE="extracted_vids.tsv"

# Write header first (optional)
zcat "$INPUT_FILE" | head -n1 | cut -f1,3,4,5,6,18,22,30 > "$OUTPUT_FILE"

# Loop through vids
while IFS= read -r vid; do
    echo "Extracting $vid..."
    zcat "$INPUT_FILE" | grep -m 1 "^$vid" | cut -f1,3,4,5,6,18,22,30 >> "$OUTPUT_FILE"
done < "$VIDS_FILE"
