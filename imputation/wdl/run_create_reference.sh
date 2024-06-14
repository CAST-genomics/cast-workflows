#!/bin/bash

region="chr15:83855424-93857434"
snp_vcf="gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz"
vntr_vcf="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_merged_samples.sorted.vcf.gz"
sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/vntr_samples_650_train.txt"
date
time python create_reference.py \
	--name vntr_reference_650 \
	--snp-vcf $snp_vcf \
	--vntr-vcf $vntr_vcf \
	--region $region \
	--samples $sample_list
date
