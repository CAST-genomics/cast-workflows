#!/bin/bash

region_file=""
samples_file=""
vcf=""
output_prefix="benchmark_1"
echo "Subsetting the target region and samples from the vcf file"
bcftools view -Oz -I -R ${regions_file} -S ${samples_file} ${vcf} > ${out_prefix}.vcf.gz && tabix -p vcf ${out_prefix}.vcf.gz
