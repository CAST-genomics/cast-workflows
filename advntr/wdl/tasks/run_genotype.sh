#!/bin/bash

# Batch size was 50 on the ACAN run
batch_size=50
mem=5
output="all_vntrs_1_sample"
region_file="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_1k_each_side.bed"
vntr_id="ALL"

time python genotype_aou.py \
  --dataset lrwgs \
  --sample-count 1 \
  --batch-size 1 \
  --output-name $output \
  --region-file $region_file \
  --mem $mem \
  --vntr-id $vntr_id
