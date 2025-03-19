#!/bin/bash

# Step 1: Group the samples in 6 differnt ancestry based on aou-predicted ancestry 
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv .
awk 'BEGIN {OFS="\t"; print "Sample_ID", "AoU_ancestry_pred"} NR>1 {print $1, $2}' ancestry_preds.tsv > flare-ga-check.txt
# awk '{print $2}' ancestry_preds.tsv | sort | uniq -c
# 56913 afr
# 45035 amr
# 5706 eas
# 133581 eur
# 942 mid
# 3217 sas

awk '$2 == "afr" {print $1}' flare-ga-check.txt > aou-afr-sample-list.txt
awk '$2 == "amr" {print $1}' flare-ga-check.txt > aou-amr-sample-list.txt
awk '$2 == "sas" {print $1}' flare-ga-check.txt > aou-sas-sample-list.txt
awk '$2 == "eas" {print $1}' flare-ga-check.txt > aou-eas-sample-list.txt
awk '$2 == "eur" {print $1}' flare-ga-check.txt > aou-eur-sample-list.txt
awk '$2 == "mid" {print $1}' flare-ga-check.txt > aou-mid-sample-list.txt
##########################################################################
##########################################################################

# Step 2: Extract and create VCF for each group of samples, splitting chunks for AFR, AMR, EUR with chink size 10K
##########################
# Download the phased VCF of all AoU samples on chr5
gsutil -u $GOOGLE_PROJECT cp gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/phased-vcf/acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz .
gsutil -u $GOOGLE_PROJECT cp gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/phased-vcf/acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz.tbi .
##########################

bcftools view -S aou-sas-sample-list.txt --force-samples \
  -Oz -o acaf_threshold_unrelated_SAS.chr5.biallelic.phased.vcf.gz \
  acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz
bcftools index --tbi acaf_threshold_unrelated_SAS.chr5.biallelic.phased.vcf.gz
echo "Creating SAS VCF complete."

##########################################################################

bcftools view -S aou-eas-sample-list.txt --force-samples \
  -Oz -o acaf_threshold_unrelated_EAS.chr5.biallelic.phased.vcf.gz \
  acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz
bcftools index --tbi acaf_threshold_unrelated_EAS.chr5.biallelic.phased.vcf.gz
echo "Creating EAS VCF complete."

##########################################################################

bcftools view -S aou-mid-sample-list.txt --force-samples \
  -Oz -o acaf_threshold_unrelated_MID.chr5.biallelic.phased.vcf.gz \
  acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz
bcftools index --tbi acaf_threshold_unrelated_MID.chr5.biallelic.phased.vcf.gz
echo "Creating MID VCF complete."

##########################################################################
##########################################################################


# Download GNU Parallel to parallel the splitting job of bcftools
#wget -qO parallel.tar.bz2 https://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2

# Extract the tarball
#tar -xf parallel.tar.bz2

#cd parallel-*

# Use GNU Parallel without installing
#/home/jupyter/workspaces/cast/parallel-20250222/src/parallel --version
# if has results: '/src/parallel' file exists, then export the path:
#export PATH=/home/jupyter/workspaces/cast/parallel-20250222/src:$PATH

#parallel --version

#echo 'export PATH=/home/jupyter/workspaces/cast/parallel-20250222/src:$PATH' >> ~/.bashrc

#source ~/.bashrc

#cd ../


##########################################################################
##########################################################################

bcftools view -S aou-eur-sample-list.txt --force-samples \
  -Oz -o acaf_threshold_unrelated_EUR.chr5.biallelic.phased.vcf.gz \
  acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz
bcftools index --tbi acaf_threshold_unrelated_EUR.chr5.biallelic.phased.vcf.gz
#####################
# Set variables 
VCF_FILE="acaf_threshold_unrelated_EUR.chr5.biallelic.phased.vcf.gz"
CHUNK_SIZE=10000
OUTPUT_DIR="ancestry-flare-EUR-run"

# Create output directory
mkdir -p $OUTPUT_DIR

# Get the total number of samples
TOTAL_SAMPLES=$(bcftools query -l $VCF_FILE | wc -l)

# Function to process each chunk
process_chunk() {
    local start=$1
    local end=$2
    local chunk_id=$3
    local vcf_file=$4
    local output_dir=$5

    # Create a temporary file for the sample IDs of this chunk
    tmp_sample_file=$(mktemp)
    bcftools query -l "$vcf_file" | sed -n "${start},${end}p" > "$tmp_sample_file"

    # Create separate VCF files for each chunk of sample IDs
    bcftools view -S "$tmp_sample_file" "$vcf_file" -Oz -o "${output_dir}/chunk_${chunk_id}.vcf.gz"

    # Clean up temporary files
    rm "$tmp_sample_file"
}

export -f process_chunk

# Generate chunk ranges and parallelize processing
seq 1 $CHUNK_SIZE $TOTAL_SAMPLES | \
    parallel -j88 --colsep ' ' '
        start={};
        end=$((start + '$CHUNK_SIZE' - 1));
        if [ $end -gt '$TOTAL_SAMPLES' ]; then end='$TOTAL_SAMPLES'; fi;
        chunk_id=$(((start - 1) / '$CHUNK_SIZE'));
        process_chunk $start $end $chunk_id '"$VCF_FILE"' '"$OUTPUT_DIR"'
    '

echo "Splitting EUR chunks complete."

##########################################################################

bcftools view -S aou-amr-sample-list.txt --force-samples \
  -Oz -o acaf_threshold_unrelated_AMR.chr5.biallelic.phased.vcf.gz \
  acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz
bcftools index --tbi acaf_threshold_unrelated_AMR.chr5.biallelic.phased.vcf.gz
#####################
# Set variables 
VCF_FILE="acaf_threshold_unrelated_AMR.chr5.biallelic.phased.vcf.gz"
CHUNK_SIZE=10000
OUTPUT_DIR="ancestry-flare-AMR-run"

# Create output directory
mkdir -p $OUTPUT_DIR

# Get the total number of samples
TOTAL_SAMPLES=$(bcftools query -l $VCF_FILE | wc -l)

# Function to process each chunk
process_chunk() {
    local start=$1
    local end=$2
    local chunk_id=$3
    local vcf_file=$4
    local output_dir=$5

    # Create a temporary file for the sample IDs of this chunk
    tmp_sample_file=$(mktemp)
    bcftools query -l "$vcf_file" | sed -n "${start},${end}p" > "$tmp_sample_file"

    # Create separate VCF files for each chunk of sample IDs
    bcftools view -S "$tmp_sample_file" "$vcf_file" -Oz -o "${output_dir}/chunk_${chunk_id}.vcf.gz"

    # Clean up temporary files
    rm "$tmp_sample_file"
}

export -f process_chunk

# Generate chunk ranges and parallelize processing
seq 1 $CHUNK_SIZE $TOTAL_SAMPLES | \
    parallel -j88 --colsep ' ' '
        start={};
        end=$((start + '$CHUNK_SIZE' - 1));
        if [ $end -gt '$TOTAL_SAMPLES' ]; then end='$TOTAL_SAMPLES'; fi;
        chunk_id=$(((start - 1) / '$CHUNK_SIZE'));
        process_chunk $start $end $chunk_id '"$VCF_FILE"' '"$OUTPUT_DIR"'
    '

echo "Splitting AMR chunks complete."

##########################################################################

bcftools view -S aou-afr-sample-list.txt --force-samples \
  -Oz -o acaf_threshold_unrelated_AFR.chr5.biallelic.phased.vcf.gz \
  acaf_threshold_unrelated.chr5.biallelic.phased.vcf.gz
bcftools index --tbi acaf_threshold_unrelated_AFR.chr5.biallelic.phased.vcf.gz
#####################
# Set variables 
VCF_FILE="acaf_threshold_unrelated_AFR.chr5.biallelic.phased.vcf.gz"
CHUNK_SIZE=10000
OUTPUT_DIR="ancestry-flare-AFR-run"

# Create output directory
mkdir -p $OUTPUT_DIR

# Get the total number of samples
TOTAL_SAMPLES=$(bcftools query -l $VCF_FILE | wc -l)

# Function to process each chunk
process_chunk() {
    local start=$1
    local end=$2
    local chunk_id=$3
    local vcf_file=$4
    local output_dir=$5

    # Create a temporary file for the sample IDs of this chunk
    tmp_sample_file=$(mktemp)
    bcftools query -l "$vcf_file" | sed -n "${start},${end}p" > "$tmp_sample_file"

    # Create separate VCF files for each chunk of sample IDs
    bcftools view -S "$tmp_sample_file" "$vcf_file" -Oz -o "${output_dir}/chunk_${chunk_id}.vcf.gz"

    # Clean up temporary files
    rm "$tmp_sample_file"
}

export -f process_chunk

# Generate chunk ranges and parallelize processing
seq 1 $CHUNK_SIZE $TOTAL_SAMPLES | \
    parallel -j88 --colsep ' ' '
        start={};
        end=$((start + '$CHUNK_SIZE' - 1));
        if [ $end -gt '$TOTAL_SAMPLES' ]; then end='$TOTAL_SAMPLES'; fi;
        chunk_id=$(((start - 1) / '$CHUNK_SIZE'));
        process_chunk $start $end $chunk_id '"$VCF_FILE"' '"$OUTPUT_DIR"'
    '

echo "Splitting AFR chunks complete."

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

wget https://raw.githubusercontent.com/browning-lab/ukb-phasing/refs/heads/master/maps/plink.chr5.GRCh38.map
gsutil -u $GOOGLE_PROJECT cp gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/reference/ucsd-ref.panel .

# Download FLARE jar
# (This example uses the official browning-lab release)
wget -O flare.jar https://faculty.washington.edu/browning/flare.jar

# Download the ref-VCF from workspace bucket
gsutil -u $GOOGLE_PROJECT cp gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/reference/hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz .
gsutil -u $GOOGLE_PROJECT cp gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/reference/hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz.tbi .

##################################################################################

# Set parameters
JAVA_OPTS="-Xmx500g"
FLARE_JAR="flare.jar"
INPUT_VCF="acaf_threshold_unrelated_SAS.chr5.biallelic.phased.vcf.gz"
REF_VCF_FILE="hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr5.GRCh38.map"
OUTPUT_PREFIX="flare_output.six-ancestry-acaf-unrelated-SAS.chr5"

# Run FLARE
java $JAVA_OPTS -jar $FLARE_JAR \
    ref=$REF_VCF_FILE \
    ref-panel=$REF_PANEL \
    gt=$INPUT_VCF \
    map=$MAP_FILE \
    probs=true \
    out=$OUTPUT_PREFIX

bcftools index --tbi flare_output.six-ancestry-acaf-unrelated-SAS.chr5.anc.vcf.gz
##############################

# run-FLARE-eas-mid.sh

# Set parameters
JAVA_OPTS="-Xmx500g"
FLARE_JAR="flare.jar"
INPUT_VCF="acaf_threshold_unrelated_EAS.chr5.biallelic.phased.vcf.gz"
REF_VCF_FILE="hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr5.GRCh38.map"
OUTPUT_PREFIX="flare_output.six-ancestry-acaf-unrelated-EAS.chr5"

# Run FLARE
java $JAVA_OPTS -jar $FLARE_JAR \
    ref=$REF_VCF_FILE \
    ref-panel=$REF_PANEL \
    gt=$INPUT_VCF \
    map=$MAP_FILE \
    probs=true \
    out=$OUTPUT_PREFIX

bcftools index --tbi flare_output.six-ancestry-acaf-unrelated-EAS.chr5.anc.vcf.gz
##############################

# Set parameters
JAVA_OPTS="-Xmx500g"
FLARE_JAR="flare.jar"
INPUT_VCF="acaf_threshold_unrelated_MID.chr5.biallelic.phased.vcf.gz"
REF_VCF_FILE="hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr5.GRCh38.map"
OUTPUT_PREFIX="flare_output.six-ancestry-acaf-unrelated-MID.chr5"

# Run FLARE
java $JAVA_OPTS -jar $FLARE_JAR \
    ref=$REF_VCF_FILE \
    ref-panel=$REF_PANEL \
    gt=$INPUT_VCF \
    map=$MAP_FILE \
    probs=true \
    out=$OUTPUT_PREFIX

bcftools index --tbi flare_output.six-ancestry-acaf-unrelated-MID.chr5.anc.vcf.gz
##################################################################################
##################################################################################

# Set variables
JAVA_OPTS="-Xmx512g"
FLARE_JAR="flare.jar"
REF_VCF_FILE="hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr5.GRCh38.map"
OUTPUT_PREFIX="ancestry-flare-EUR-run/flare_output.six-ancestry-acaf-unrelated-EUR.chr5"
OUTPUT_DIR="ancestry-flare-EUR-run"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to run FLARE on each chunk
process_chunk_with_flare() {
    # Process each chunk in numerical order
    find ${OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | while read chunk_file; do
        chunk_id=$(basename $chunk_file | sed 's/chunk_\([0-9]*\).vcf.gz/\1/')
        
        # Define log file for each chunk
        LOG_FILE="${OUTPUT_PREFIX}.chunk_${chunk_id}.log"
        
        # Run FLARE for the chunk and capture log
        java $JAVA_OPTS -jar $FLARE_JAR \
            ref=$REF_VCF_FILE \
            ref-panel=$REF_PANEL \
            gt=$chunk_file \
            map=$MAP_FILE \
            probs=true \
            em=true \
            out=${OUTPUT_PREFIX}.chunk_${chunk_id} &> "$LOG_FILE"
    done

    echo "FLARE EUR processing complete."
}

# Function to index all .anc.vcf.gz files in the output directory
index_anc_vcf_files() {
    echo "Indexing all .anc.vcf.gz files..."
    find $OUTPUT_DIR -name "*.anc.vcf.gz" | while read vcf_file; do
        bcftools index --tbi "$vcf_file"
    done
    echo "Indexing EUR complete."
}

# Run FLARE processing
process_chunk_with_flare

# Run indexing for all .anc.vcf.gz files
index_anc_vcf_files

##########################################################################

# Set variables
JAVA_OPTS="-Xmx512g"
FLARE_JAR="flare.jar"
REF_VCF_FILE="hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr5.GRCh38.map"
OUTPUT_PREFIX="ancestry-flare-AFR-run/flare_output.six-ancestry-acaf-unrelated-AFR.chr5"
OUTPUT_DIR="ancestry-flare-AFR-run"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to run FLARE on each chunk
process_chunk_with_flare() {
    # Process each chunk in numerical order
    find ${OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | while read chunk_file; do
        chunk_id=$(basename $chunk_file | sed 's/chunk_\([0-9]*\).vcf.gz/\1/')
        
        # Define log file for each chunk
        LOG_FILE="${OUTPUT_PREFIX}.chunk_${chunk_id}.log"
        
        # Run FLARE for the chunk and capture log
        java $JAVA_OPTS -jar $FLARE_JAR \
            ref=$REF_VCF_FILE \
            ref-panel=$REF_PANEL \
            gt=$chunk_file \
            map=$MAP_FILE \
            probs=true \
            em=true \
            out=${OUTPUT_PREFIX}.chunk_${chunk_id} &> "$LOG_FILE"
    done

    echo "FLARE AFR processing complete."
}

# Function to index all .anc.vcf.gz files in the output directory
index_anc_vcf_files() {
    echo "Indexing all .anc.vcf.gz files..."
    find $OUTPUT_DIR -name "*.anc.vcf.gz" | while read vcf_file; do
        bcftools index --tbi "$vcf_file"
    done
    echo "Indexing AFR complete."
}

# Run FLARE processing
process_chunk_with_flare

# Run indexing for all .anc.vcf.gz files
index_anc_vcf_files


##########################################################################

# Set variables
JAVA_OPTS="-Xmx512g"
FLARE_JAR="flare.jar"
REF_VCF_FILE="hgdp-1kgp-ref_chr5_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr5.GRCh38.map"
OUTPUT_PREFIX="ancestry-flare-AMR-run/flare_output.six-ancestry-acaf-unrelated-AMR.chr5"
OUTPUT_DIR="ancestry-flare-AMR-run"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to run FLARE on each chunk
process_chunk_with_flare() {
    # Process each chunk in numerical order
    find ${OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | while read chunk_file; do
        chunk_id=$(basename $chunk_file | sed 's/chunk_\([0-9]*\).vcf.gz/\1/')
        
        # Define log file for each chunk
        LOG_FILE="${OUTPUT_PREFIX}.chunk_${chunk_id}.log"
        
        # Run FLARE for the chunk and capture log
        java $JAVA_OPTS -jar $FLARE_JAR \
            ref=$REF_VCF_FILE \
            ref-panel=$REF_PANEL \
            gt=$chunk_file \
            map=$MAP_FILE \
            probs=true \
            em=true \
            out=${OUTPUT_PREFIX}.chunk_${chunk_id} &> "$LOG_FILE"
    done

    echo "FLARE AMR processing complete."
}

# Function to index all .anc.vcf.gz files in the output directory
index_anc_vcf_files() {
    echo "Indexing all .anc.vcf.gz files..."
    find $OUTPUT_DIR -name "*.anc.vcf.gz" | while read vcf_file; do
        bcftools index --tbi "$vcf_file"
    done
    echo "Indexing AMR complete."
}

# Run FLARE processing
process_chunk_with_flare

# Run indexing for all .anc.vcf.gz files
index_anc_vcf_files

##################################################################################
##################################################################################
##################################################################################
##################################################################################

# Define output file
OUTPUT_GLOBAL="flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.global.anc.gz"

# Find all global ancestry files to merge
FILES_TO_MERGE=$(find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run -name "*.global.anc.gz")
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-EAS.chr5.global.anc.gz"
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-MID.chr5.global.anc.gz"
FILES_TO_MERGE+=" flare_output.six-ancestry-acaf-unrelated-SAS.chr5.global.anc.gz"

# Check if all files exist
echo "Merging the following files:"
echo "$FILES_TO_MERGE"

# Extract header from the first file and merge all while removing duplicate headers
(zcat flare_output.six-ancestry-acaf-unrelated-EAS.chr5.global.anc.gz | head -n 1;
 for file in $FILES_TO_MERGE; do
     zcat "$file" | tail -n +2
 done) | gzip > "$OUTPUT_GLOBAL"

zcat flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.global.anc.gz | \
awk 'NR==1 {print $1, $7, $3, $5, $4, $2, $6; next} {print $1, $7, $3, $5, $4, $2, $6}' OFS="\t" | \
bgzip > flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.global.anc.tmp.gz

# Replace the original file safely
mv flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.global.anc.tmp.gz \
   flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.global.anc.gz

echo "Merged global ancestry file created: $OUTPUT_GLOBAL"

##################################################################################

# Define output file names
OUTPUT_MODEL="flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.model"
OUTPUT_LOG="flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.log"

# Find all .model files and merge them
echo "Merging .model files..."
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run -name "*.model" | sort -V > model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr5.model" >> model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr5.model" >> model_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr5.model" >> model_files.txt

# Concatenate all model files into one
cat $(cat model_files.txt) > "$OUTPUT_MODEL"
rm model_files.txt

echo "Merged .model file created: $OUTPUT_MODEL"

# Find all .log files and merge them
echo "Merging .log files..."
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run -name "*.log" | sort -V > log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr5.log" >> log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr5.log" >> log_files.txt
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr5.log" >> log_files.txt

# Concatenate all log files into one
cat $(cat log_files.txt) > "$OUTPUT_LOG"
rm log_files.txt

echo "Merged .log file created: $OUTPUT_LOG"

##################################################################################

# Define output file
OUTPUT_VCF="flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.anc.vcf.gz"

# Create a temporary file list
FILE_LIST="vcf_file_list.txt"

# Find all VCF files and write them to a file list
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run -name "*.anc.vcf.gz" | sort -V > "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr5.anc.vcf.gz" >> "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr5.anc.vcf.gz" >> "$FILE_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr5.anc.vcf.gz" >> "$FILE_LIST"

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

bcftools index --tbi flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.anc.vcf.gz

##################################################################################
##################################################################################
##################################################################################
##################################################################################

# Upload files to workspace bucket:
gsutil -u $GOOGLE_PROJECT cp flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.global.anc.gz gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/
gsutil -u $GOOGLE_PROJECT cp flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.model gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/
gsutil -u $GOOGLE_PROJECT cp flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.log gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/
gsutil -u $GOOGLE_PROJECT cp flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.anc.vcf.gz gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/
gsutil -u $GOOGLE_PROJECT cp flare_output.seperated-ancestry-run-acaf-unrelated.chr5.merged.anc.vcf.gz.tbi gs://fc-secure-5efb5f8d-7703-495b-8dc7-3d90605152e5/ancestry-calling/


