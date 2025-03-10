# wgonzalezrivera
# last update: March 10, 2025
# In this file, I will write down how flare segments are being retrieve. We will have chunk_9 as test.

# STEP 1: Get .tsv similar to .msp file in gnomix/rfmix
# Example: chrom\t pos\t samp(x+1)-hap1\t samp(x+1)-hap2 #TODO: ADD SNP_ID

#!/bin/bash

# Define input and output files
VCF="/home/jupyter/flare_chr_segments_weighted_ancestry_proportions/flare_output.6-ancestry-acaf-unrelated-all.chr21.chunk_9.anc.vcf.gz"
OUTPUT="chr21_chunk_9_testing_msp.tsv"

# Step 1: Generate the header row
{
  # Print the fixed header columns
  echo -ne "chr\tpos\t"
  
  # Extract sample names and append "-h1" and "-h2" for each sample
  bcftools query -l "$VCF" | awk '{printf "%s-h1\t%s-h2\t", $1, $1}'
  
  # Add a newline to complete the header
  echo ""
} > "$OUTPUT"

# Step 2: Extract data, filter rows with missing ancestry data, and sort by POS
bcftools query -f "%CHROM\t%POS[\t%AN1\t%AN2]\n" "$VCF" |
  awk '!/\t\.\t|\t\.$/' |  # Remove rows with missing ancestry data (denoted by ".")
  sort -k2,2n -S 32G --parallel=8 >> "$OUTPUT"  # Sort by the second column (POS) numerically

# Testing if this code works... 
