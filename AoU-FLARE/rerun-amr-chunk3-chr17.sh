#!/bin/bash

# Special case handling: rerun FLARE local ancestry calling on Chr17 AMR-chunk_3 & Chr18 AMR-chunk_3
# Standard workflow failed on Chr17 AMR-chunk_3 (chunk sample size =10K) & Chr18 AMR-chunk_3 (chunk sample size =10K)
# Issue solved: merge AMR-chunk_3 and AMR-chunk_4 (chunk sample size =11,018), then rerun FLARE on this merged AMR_chunk_3_4 VCF


# Step 1: Merge AMR-chunk_3 and AMR-chunk_4

# Index the input VCFs if not already done
tabix -p vcf flare-AMR-run/chunk_3.vcf.gz
tabix -p vcf flare-AMR-run/chunk_4.vcf.gz

# Merge using bcftools merge
bcftools merge flare-AMR-run/chunk_3.vcf.gz flare-AMR-run/chunk_4.vcf.gz -Oz -o flare-AMR-run/merged_chunk_3_4.vcf.gz

# Index the merged VCF
tabix -p flare-AMR-run/vcf merged_chunk_3_4.vcf.gz


###################################################################################################################
###################################################################################################################

# Step 2: Rerun FLARE on the merged AMR-chunk_3_4 VCF file
# Set parameters
JAVA_OPTS="-Xmx500g"
FLARE_JAR="flare.jar"
INPUT_VCF="ancestry-flare-AMR-run/merged_chunk_3_4.vcf.gz"
REF_VCF_FILE="hgdp-1kgp-ref_chr17_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr17.GRCh38.map"
OUTPUT_PREFIX="ancestry-flare-AMR-run/flare_output.six-ancestry-acaf-unrelated-AMR.chr17.chunk_3_4"

# Run FLARE
java $JAVA_OPTS -jar $FLARE_JAR \
    ref=$REF_VCF_FILE \
    ref-panel=$REF_PANEL \
    gt=$INPUT_VCF \
    map=$MAP_FILE \
    probs=true \
    em=true \
    out=$OUTPUT_PREFIX

bcftools index --tbi ancestry-flare-AMR-run/flare_output.six-ancestry-acaf-unrelated-AMR.chr17.chunk_3_4.anc.vcf.gz


###################################################################################################################
###################################################################################################################

# Step 3: List and Merge all chunk output files
# move the previous AMR-chunk_3 FLARE outputs to the folder 'pre-run'
# mkdir -p pre-run


# Define output file
OUTPUT_GLOBAL="flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.global.anc.gz"

# Find all global ancestry files to merge
FILES_TO_MERGE=$(find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.global.anc.gz" ! -path "*/pre-run/*")
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-EAS.chr17.global.anc.gz"
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-MID.chr17.global.anc.gz"
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-SAS.chr17.global.anc.gz"

# Check if all files exist
echo "Merging the following files:"
echo "$FILES_TO_MERGE"

# Extract header from the first file and merge all while removing duplicate headers
(zcat flare_output.six-ancestry-acaf-unrelated-EAS.chr17.global.anc.gz | head -n 1;
 for file in $FILES_TO_MERGE; do
     zcat "$file" | tail -n +2
 done) | gzip > "$OUTPUT_GLOBAL"

zcat flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.global.anc.gz | \
awk 'NR==1 {print $1, $7, $3, $5, $4, $2, $6; next} {print $1, $7, $3, $5, $4, $2, $6}' OFS="\t" | \
bgzip > flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.global.anc.tmp.gz

# Replace the original file safely
mv flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.global.anc.tmp.gz \
   flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.global.anc.gz

echo "Merged global ancestry file created: $OUTPUT_GLOBAL"

##################################################################################

# Define output file names
OUTPUT_MODEL="flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.model"
OUTPUT_LOG="flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.log"

# Find all .model files and merge them
echo "Merging .model files..."
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.model" ! -path "*/pre-run/*" | sort -V > model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr17.model" >> model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr17.model" >> model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr17.model" >> model_files.txt

# Concatenate all model files into one
cat $(cat model_files.txt) > "$OUTPUT_MODEL"
rm model_files.txt

echo "Merged .model file created: $OUTPUT_MODEL"



# Find all .log files (excluding any inside pre-run subdirectories)
echo "Merging .log files..."
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.log" ! -path "*/pre-run/*" | sort -V > log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr17.log" >> log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr17.log" >> log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr17.log" >> log_files.txt

# Concatenate all log files into one
cat $(cat log_files.txt) > "$OUTPUT_LOG"
rm log_files.txt

echo "Merged .log file created: $OUTPUT_LOG"


##################################################################################

# Define output file
OUTPUT_VCF="flare_output.seperated-ancestry-run-acaf-unrelated.chr17.merged.anc.vcf.gz"

# Create a temporary file list
FILE_LIST="vcf_file_list.txt"

# Find all VCF files and write them to a file list
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.anc.vcf.gz" ! -path "*/pre-run/*" | sort -V > "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr17.anc.vcf.gz" >> "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr17.anc.vcf.gz" >> "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr17.anc.vcf.gz" >> "$FILE_LIST"

# Display the files to be merged
echo "Merging the following VCF files:"
cat "$FILE_LIST"

# Run bcftools merge
bcftools merge --file-list "$FILE_LIST" -Oz -o "$OUTPUT_VCF"

# Index the merged VCF
bcftools index --tbi "$OUTPUT_VCF"

# Clean up temporary file
rm "$FILE_LIST"

echo "Merged VCF file created: $OUTPUT_VCF"
