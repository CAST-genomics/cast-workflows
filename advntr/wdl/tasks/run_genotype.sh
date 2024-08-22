#!/bin/bash

# Batch size was 50 on the ACAN run
batch_size=50
mem=20
output="all_vntrs_1_sample"
region_file="$WORKSPACE_BUCKET/saraj/p_vntrs_g_vntrs/p_vntrs_g_vntrs_merged.bed"
vntr_id="ALL"

time python genotype_aou.py \
  --dataset lrwgs \
  --sample-count 4 \
  --batch-size 2 \
  --output-name $output \
  --region-file $region_file \
  --mem $mem \
  --vntr-id $vntr_id
