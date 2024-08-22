#!/bin/bash

chr="chr15"

#mem=20
mem=4
#window=20
window=3
#regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_50m.bed"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_10m.bed"

snp_vcf="gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz"
#vntr_vcf="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_merged_samples.sorted.vcf.gz"
vntr_vcf="$WORKSPACE_BUCKET/saraj/acan_data_test/merged_4_samples_batches.sorted.vcf.gz"
#sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/vntr_samples_650_train.txt"
#sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/vntr_samples_20_small_test.txt"
sample_list="$WORKSPACE_BUCKET/saraj/acan_data_test/vntr_samples_4_small_test.txt"
map="$WORKSPACE_BUCKET/saraj/genetic_map/plink.${chr}_w_chr.GRCh38.map"

date
time python create_reference.py \
	--name vntr_ref_test \
	--snp-vcf $snp_vcf \
	--vntr-vcf $vntr_vcf \
	--regions $regions \
	--samples $sample_list \
        --window $window \
	--chrom $chr \
	--mem $mem \
	--cromwell \
	--map $map

date
