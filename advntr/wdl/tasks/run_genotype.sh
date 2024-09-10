#!/bin/bash

# Batch size was 50 on the ACAN run
batch_size=300
#batch_size=2
mem=60
batch_idx=0
output="p_g_vntrs_batch_all_repeat_missing_id"
#output="p_vntrs_merge_id_repeat_2"
#output="p_g_vntrs_gte_500bp_test"
region_file="$WORKSPACE_BUCKET/saraj/p_vntrs_g_vntrs/p_vntrs_g_vntrs_merged.bed"
#region_file="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_1k_each_side.bed"
#vntr_id="ALL"
#vntr_id="290964"
# p_vntrs without the 6 p_vntrs with the longest runtime (> 100s)
# IDs with long runtime: 983809, 165038, 731268, 385941, 666245, 123860
p_vntr_ids_fast="290964,4237,674126,628764,376172,450626,49244,239917,705182,983810,983811,75781,983812,983813,826396,331665,731278,583186,493003,122595,886134,45940,181577,666226,666244,915594,492236,983814,983815,527656,983816,339765,983817,358459,936688,936686,944729,983820,983821,746692,983808,665101"
p_vntr_ids="290964,4237,674126,628764,376172,450626,49244,239917,705182,983810,983811,75781,983812,983813,826396,331665,731278,583186,493003,122595,886134,45940,181577,666226,666244,915594,492236,983814,983815,527656,983816,339765,983817,358459,936688,936686,944729,983820,983821,746692,983808,665101,983809,165038,731268,385941,666245,123860"
vntr_id="ALL"

time python genotype_aou.py \
  --dataset lrwgs \
  --sample-count 1027 \
  --batch-size $batch_size \
  --output-name $output \
  --region-file $region_file \
  --mem $mem \
  --only-batch-index $batch_idx \
  --vntr-id $vntr_id
