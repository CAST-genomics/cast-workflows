#!/bin/bash

vcf_file=$1
if [ ! -e ${vcf_file}.tbi ]; then
  tabix -p vcf $vcf_file
fi
bcftools view -i 'ID="."' $vcf_file | grep -v "^#" | wc -l
