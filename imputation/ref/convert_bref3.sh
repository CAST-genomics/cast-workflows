#!/bin/bash

#set up gc bucket storage"
gcs="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/"
chr=$1

#dowaload bref3.jar
bref="bref3.27May24.118.jar"
# Check if the file exists locally
if [ ! -f "$bref" ]; then
    echo "File does not exist. Downloading..."
    wget https://faculty.washington.edu/browning/beagle/bref3.27May24.118.jar 
else
    echo "File already exists. Continuing..."
fi

#convert vcf to bref3
echo "converting vcf to bref3 format"
zcat ${chr}_final_SNP_merged_additional_TRs.vcf.gz | java -jar $bref > ${chr}_final_SNP_merged_additional_TRs.bref3

echo "Done converting"

gsutil cp ${chr}_final_SNP_merged_additional_TRs.bref3 $gcs
echo "Successfully transfer bref to google cloud bucket"
 



