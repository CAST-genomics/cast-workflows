#!/bin/bash

## We use PheTK to run phewas. For more information on installation and documentation see https://github.com/nhgritctran/PheTK/tree/main?tab=readme-ov-file

## To install run:
#pip install PheTK --upgrade

#./run_aou_pheTK.sh ../targetTR/strsets/ukb_finemapped_hg38/ukb_finemapped_hg38_batch2.bed ~/compare_hipstr_imputation/UKB_finemapped-batch2.filtered.sorted.vcf.gz

# Path to the region file and tr_vcf file
region_file="$1"
tr_vcf="$2"

# Ensure the region file exists
if [[ ! -f "$region_file" ]]; then
    echo "Region file does not exist."
    exit 1
fi

# Iterate over each line in the region file
while IFS= read -r line; do
    # Remove leading and trailing whitespaces
    line=$(echo "$line" | xargs)
    region=$(echo "$line" | awk '{print $1 ":" $2 "-" $3}')
    
    # Run the Python script with the provided arguments
    python run_pheTK.py --region "$region" --tr-vcf "$tr_vcf" --n-threads 12 --out "$region_pheTK"  
done < "$region_file"
