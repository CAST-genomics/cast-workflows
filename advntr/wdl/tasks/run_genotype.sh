#!/bin/bash

# Batch size was 50 on the ACAN run
batch_size=50
mem=20
output="all_vntrs_10_samples"
region_file="$WORKSPACE_BUCKET/saraj/p_vntrs_g_vntrs/p_vntrs_g_vntrs_merged.bed"
#region_file="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_1k_each_side.bed"
vntr_id="ALL"
#vntr_id="290964"

time python genotype_aou.py \
  --dataset lrwgs \
  --sample-count 10 \
  --batch-size 5 \
  --output-name $output \
  --region-file $region_file \
  --mem $mem \
  --vntr-id $vntr_id
