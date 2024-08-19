#!/bin/bash

chr="chr15"

#window=20
window=2
#regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_50m.bed"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_4m.bed"

snp_vcf="gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz"
vntr_vcf="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_merged_samples.sorted.vcf.gz"
sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/vntr_samples_650_train.txt"
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
	--mem 20 \
	--cromwell \
	--map $map

date
