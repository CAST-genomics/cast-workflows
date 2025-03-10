# wgonzalezrivera
# last update: March 10, 2025
# In this file, I will write down how flare segments are being retrieve. We will have chunk_9 as test.

# STEP 1: Get .tsv similar to .msp file in gnomix/rfmix
# Example: chrom\t pos\t samp(x+1)-hap1\t samp(x+1)-hap2 #TODO: ADD SNP_ID

#!/bin/bash

# Define input and output files
VCF="path_to_file/flare_output.6-ancestry-acaf-unrelated-all.chr21.chunk_9.anc.vcf.gz"
OUTPUT="chr21_chunk_9_msp.tsv"

# Step 1: Generate the header row
{
  # Print the fixed header columns
  echo -ne "chr\tpos\t"
  
  # Extract sample names and append "-h1" and "-h2" for each sample
  bcftools query -l "$VCF" | awk '{printf "%s-h1\t%s-h2\t", $1, $1}'
  
  # Add a newline to complete the header
  echo ""
} > "$OUTPUT"

# Step 2: Extract data, filter rows with missing or invalid ancestry data, and sort by POS
bcftools query -f "%CHROM\t%POS[\t%AN1\t%AN2]\n" "$VCF" |
  awk '
    # Remove rows with missing ancestry data (denoted by ".")
    !/\t\.\t|\t\.$/ {
      # Check if columns 3 and beyond are valid integers
      valid = 1
      for (i = 3; i <= NF; i++) {
        if ($i !~ /^-?[0-9]+$/) {
          valid = 0
          break
        }
      }
      # Print the row only if all columns are valid
      if (valid) {
        print
      }
    }
  ' |
  sort -k2,2n -S 32G --parallel=8 >> "$OUTPUT"  # Sort by the second column (POS) numerically

# Step 3: running process_segments.py to get segments\t labels\t lengths\t 
# Will try process_segments.py before uploading to github

./process_segments.py chr21_chunk_9_testing_msp.tsv chr21_chunk_9_testing_segments.tsv

# Step 4: getting the ancestry proportions by weights

import pandas as pd
import numpy as np
from collections import defaultdict

# Load DataFrame
df = pd.read_csv("/home/jupyter/flare_chr_segments_weighted_ancestry_proportions/chr21_chunk_9_segments.tsv", sep="\t")

# Convert "Lengths" column to lists of integers
df["Lengths"] = df["Lengths"].apply(lambda x: list(map(int, x.split(","))) if isinstance(x, str) else [int(x)])

# Initialize dictionary
samp_to_segments = defaultdict(lambda: [[], []])

# Populate dictionary
for _, row in df.iterrows():
    sample = row["Sample"]
    labels = row["Labels"].split(",") if isinstance(row["Labels"], str) else [row["Labels"]]
    lengths = row["Lengths"]

    samp_to_segments[sample][0].extend(labels)   # Append labels as list
    samp_to_segments[sample][1].extend(lengths)  # Append lengths as list

# Mapping of numeric ancestry codes to population labels
num_to_pop = {0: "MIDDLE_EAST", 1: "AMERICA", 2: "EUR", 3: "EAS", 4: "SAS", 5: "AFR"} #THIS CAN CHANGE, NO PARSE

# Prepare output file
output_file = "/home/jupyter/flare_chr_segments_weighted_ancestry_proportions/flare-summstats_chunk_9_chr21.tab"
with open(output_file, "w") as f:
    # Prepare the header
    header = ["SAMPLE"]
    for num in num_to_pop:
        pop = num_to_pop[num]
        for metric in ["perc", "len_mean", "len_median", "num_segments"]:
            header.append(f"{pop}-{metric}")
    f.write("\t".join(header) + "\n")

    # Compute stats for each sample
    for sample in samp_to_segments.keys():
        outitems = [sample]
        labels = samp_to_segments[sample][0]
        lengths = samp_to_segments[sample][1]
        total_len = sum(lengths)
        
        for num in num_to_pop:
            popidx = [i for i in range(len(labels)) if labels[i] == str(num)]  # Ensure labels are strings
            if len(popidx) == 0:
                outitems.extend([0, "NA", "NA", 0])
                continue
            
            popsegs = [lengths[i] for i in popidx]
            perc = sum(popsegs) / total_len
            len_mean = np.mean(popsegs)
            len_median = np.median(popsegs)
            num_segments = len(popsegs)
            
            outitems.extend([perc, len_mean, len_median, num_segments])
        
        f.write("\t".join(map(str, outitems)) + "\n")

print(f"File saved: {output_file}")
