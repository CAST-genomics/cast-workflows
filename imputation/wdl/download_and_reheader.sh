#!/bin/bash

chr=$1
outfile_tmp="imputed_${chr}.sorted.annotated.vcf.gz"
outfile="imputed_${chr}.annotated.rh.vcf.gz"

# Download from the bucket
gsutil cp $WORKSPACE_BUCKET/saraj/imputation_output/${chr}/imputed_${chr}.sorted.annotated.vcf.gz* .

tabix -p vcf ${outfile_tmp}
bcftools view -h ${outfile_tmp} | head -n 4 > header.txt
echo "##source=adVNTR ver. 1.5.0" >> header.txt
bcftools view -h ${outfile_tmp} | tail -n +5 >> header.txt
bcftools reheader -h header.txt ${outfile_tmp} | bcftools view -Oz > ${outfile}
tabix -p vcf ${outfile}
