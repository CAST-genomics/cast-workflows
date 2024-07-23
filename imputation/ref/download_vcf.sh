#!/bin/bash

#./download_vcf.sh chr1
# Define the base URL where files are located
url="https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/"
gcs="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/"
chr=$1



# Construct the URL for the current chromosome
vcf="${chr}_final_SNP_merged_additional_TRs.vcf.gz"
vcf_index="${chr}_final_SNP_merged_additional_TRs.vcf.gz.tbi"
    
# Download the file using wget
echo "Downloading ${vcf} ..."
wget "${url}${vcf}"
wget "${url}${vcf_index}"

echo "Move ref vcf files to google cloud"
gsutil cp $vcf $gcs
gsutil cp $vcf_index $gcs

echo "Download complete."

