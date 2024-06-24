#!/bin/bash

# for 70% train and 30% test:
#num_train=650
#num_test=275

# for 80% train and 20% test:
num_train=740
num_test=185
# NOTE: each time running shuf, there is a new shuffle, therefore, new set of samples. So don't run it unless necessary,
#cat vntr_samples_clean_sorted.txt | shuf | head -n $num_train > vntr_samples_${num_train}_train.txt
#cat vntr_samples_clean_sorted.txt | grep -v -f vntr_samples_${num_train}_train.txt > vntr_samples_${num_test}_test.txt

# test the samples are split correctly
cat vntr_samples_${num_train}_train.txt  vntr_samples_${num_test}_test.txt | sort --uniq > test.txt

echo "Num lines in the merged file (should be 925): $(wc -l test.txt)"
echo "Num lines in train (should be $num_train): $(wc -l vntr_samples_${num_train}_train.txt)"
echo "Num lines in test (should be $num_test): $(wc -l vntr_samples_${num_test}_test.txt)"
echo "difference in original and reconstructed sample size (should be 0) : $(comm -3 test.txt vntr_samples_clean_sorted.txt | wc -l)"

# If a tag e.g. END is missing in the header, run:
#bcftools view -h ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz > header.txt
# change the header.txt, i.e. get the ENG tag from merged_samples.sorted.vcf.gz.
# line to add for the END tag:
# ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">
# ##INFO=<ID=VID,Number=1,Type=Integer,Description="VNTR ID">
# ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat motif">
# ##INFO=<ID=RC,Number=1,Type=Integer,Description="Reference repeat unit count">
#bcftools reheader -h -p ref_phased_reheader_chr15_acan_vntr.vcf.gz ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz

tabix -f -p vcf ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz

bcftools view -Oz -S vntr_samples_${num_train}_train.txt \
	ref_phased_reheader_chr15_acan_vntr.vcf.gz \
	> phased_ACAN_vntr_snp_${num_train}_samples.sorted.vcf.gz && \
tabix -p vcf phased_ACAN_vntr_snp_${num_train}_samples.sorted.vcf.gz
