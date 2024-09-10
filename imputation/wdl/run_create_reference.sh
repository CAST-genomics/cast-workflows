#!/bin/bash

#for chr_idx in $(seq 1 22); do
for chr_idx in "22" ; do
  chr="chr${chr_idx}"
  echo "Running for $chr"
  mem=60
  #mem=4
  #window=20
  window=10
  #window=3
  #regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_50m.bed"
  #regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_10m.bed"
  regions="$WORKSPACE_BUCKET/saraj/p_vntrs_g_vntrs/p_vntrs_g_vntrs_50mbp_merged.bed"
  # To add flanking regions to the bed file
  #bedtools slop -i p_vntrs_g_vntrs_merged.bed -g hg38.chrom.sizes -b 50000000 > p_vntrs_g_vntrs_50mbp.bed
  # bedtools merge -i p_vntrs_g_vntrs_50mbp.bed > p_vntrs_g_vntrs_50mbp_merged.bed
  snp_vcf="gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz"
  #vntr_vcf="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_merged_samples.sorted.vcf.gz"
  #vntr_vcf="$WORKSPACE_BUCKET/saraj/acan_data_test/merged_4_samples_batches.sorted.vcf.gz"
  vntr_vcf="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/p_g_vntrs/merged_all_lrwgs_p_g_vntrs.sorted.vcf.gz"
  #sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/vntr_samples_650_train.txt"
  #sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/vntr_samples_20_small_test.txt"
  sample_list="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/p_g_vntrs/sample_ids_lrwgs.sorted.txt"
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
  	--map $map
  
  date
done
