#!/bin/bash

# Special case handling: rerun FLARE local ancestry calling on Chr3-EUR-chunk_2
# Standard workflow failed on Chr3-EUR-chunk_2 (chunk sample size =10K)
# Issue solved: Split Chr3-EUR-chunk_2 into 2 chunks: 
# Chr3-EUR-chunk_2.part1 (chunk sample size =5K) & Chr3-EUR-chunk_2.part2 (chunk sample size =5K), then rerun FLARE on each VCF


# Step 1: Extract all sample names
bcftools query -l chunk_2.vcf.gz > chunk_2.samples.txt

# Step 2: Split the sample list into two equal parts
head -n 5000 chunk_2.samples.txt > chunk_2.part1.samples.txt
tail -n +5001 chunk_2.samples.txt > chunk_2.part2.samples.txt

# Step 3: Extract the first 5000 samples into a new VCF
bcftools view -S chunk_2.part1.samples.txt -Oz -o chunk_2.part1.vcf.gz chunk_2.vcf.gz
tabix -p vcf chunk_2.part1.vcf.gz

# Step 4: Extract the second 5000 samples into another new VCF
bcftools view -S chunk_2.part2.samples.txt -Oz -o chunk_2.part2.vcf.gz chunk_2.vcf.gz
tabix -p vcf chunk_2.part2.vcf.gz

# Optional: Verify sample counts
echo "Samples in part 1:"
bcftools query -l chunk_2.part1.vcf.gz | wc -l

echo "Samples in part 2:"
bcftools query -l chunk_2.part2.vcf.gz | wc -l


###################################################################################################################
###################################################################################################################
# Set parameters
JAVA_OPTS="-Xmx500g"
FLARE_JAR="flare.jar"
INPUT_VCF="ancestry-flare-EUR-run/chunk_2.part1.vcf.gz"
REF_VCF_FILE="hgdp-1kgp-ref_chr3_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr3.GRCh38.map"
OUTPUT_PREFIX="ancestry-flare-EUR-run/flare_output.six-ancestry-acaf-unrelated-EUR.chr3.chunk_2.part1"

# Run FLARE
java $JAVA_OPTS -jar $FLARE_JAR \
    ref=$REF_VCF_FILE \
    ref-panel=$REF_PANEL \
    gt=$INPUT_VCF \
    map=$MAP_FILE \
    probs=true \
    em=true \
    out=$OUTPUT_PREFIX

bcftools index --tbi ancestry-flare-EUR-run/flare_output.six-ancestry-acaf-unrelated-EUR.chr3.chunk_2.part1.anc.vcf.gz


# Set parameters
JAVA_OPTS="-Xmx500g"
FLARE_JAR="flare.jar"
INPUT_VCF="ancestry-flare-EUR-run/chunk_2.part2.vcf.gz"
REF_VCF_FILE="hgdp-1kgp-ref_chr3_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr3.GRCh38.map"
OUTPUT_PREFIX="ancestry-flare-EUR-run/flare_output.six-ancestry-acaf-unrelated-EUR.chr3.chunk_2.part2"

# Run FLARE
java $JAVA_OPTS -jar $FLARE_JAR \
    ref=$REF_VCF_FILE \
    ref-panel=$REF_PANEL \
    gt=$INPUT_VCF \
    map=$MAP_FILE \
    probs=true \
    em=true \
    out=$OUTPUT_PREFIX

bcftools index --tbi ancestry-flare-EUR-run/flare_output.six-ancestry-acaf-unrelated-EUR.chr3.chunk_2.part2.anc.vcf.gz


###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

# Define output file
OUTPUT_GLOBAL="flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.global.anc.gz"

# Find all global ancestry files to merge
FILES_TO_MERGE=$(find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.global.anc.gz" ! -path "*/pre-run/*")
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-EAS.chr3.global.anc.gz"
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-MID.chr3.global.anc.gz"
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-SAS.chr3.global.anc.gz"

# Check if all files exist
echo "Merging the following files:"
echo "$FILES_TO_MERGE"

# Extract header from the first file and merge all while removing duplicate headers
(zcat flare_output.six-ancestry-acaf-unrelated-EAS.chr3.global.anc.gz | head -n 1;
 for file in $FILES_TO_MERGE; do
     zcat "$file" | tail -n +2
 done) | gzip > "$OUTPUT_GLOBAL"

zcat flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.global.anc.gz | \
awk 'NR==1 {print $1, $7, $3, $5, $4, $2, $6; next} {print $1, $7, $3, $5, $4, $2, $6}' OFS="\t" | \
bgzip > flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.global.anc.tmp.gz

# Replace the original file safely
mv flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.global.anc.tmp.gz \
   flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.global.anc.gz

echo "Merged global ancestry file created: $OUTPUT_GLOBAL"

##################################################################################

# Define output file names
OUTPUT_MODEL="flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.model"
OUTPUT_LOG="flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.log"

# Find all .model files and merge them
echo "Merging .model files..."
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.model" ! -path "*/pre-run/*" | sort -V > model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr3.model" >> model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr3.model" >> model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr3.model" >> model_files.txt

# Concatenate all model files into one
cat $(cat model_files.txt) > "$OUTPUT_MODEL"
rm model_files.txt

echo "Merged .model file created: $OUTPUT_MODEL"



# Find all .log files (excluding any inside pre-run subdirectories)
echo "Merging .log files..."
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.log" ! -path "*/pre-run/*" | sort -V > log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr3.log" >> log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr3.log" >> log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr3.log" >> log_files.txt

# Concatenate all log files into one
cat $(cat log_files.txt) > "$OUTPUT_LOG"
rm log_files.txt

echo "Merged .log file created: $OUTPUT_LOG"


##################################################################################

# Define output file
OUTPUT_VCF="flare_output.seperated-ancestry-run-acaf-unrelated.chr3.merged.anc.vcf.gz"

# Create a temporary file list
FILE_LIST="vcf_file_list.txt"

# Find all VCF files and write them to a file list
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run \
    -type f -name "*.anc.vcf.gz" ! -path "*/pre-run/*" | sort -V > "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr3.anc.vcf.gz" >> "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr3.anc.vcf.gz" >> "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr3.anc.vcf.gz" >> "$FILE_LIST"

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
