#!/bin/bash

vcf_file=$1
tabix -p vcf $vcf_file
bcftools view -i 'ID="."' $vcf_file | grep -v "^#" | wc -l
