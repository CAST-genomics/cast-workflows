#!/bin/bash

chr="chr2"

# PVNTR ref
ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/p_g_vntrs/phased/${chr}/vntr_ref_${chr}.sorted.annotated.vcf.gz"


## Run with AoU data
# VCF file with all samples
#gt="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.${chr}.vcf.bgz"
gt="$WORKSPACE_BUCKET/acaf_batches/${chr}"

# VCF file for all aou samples, only acan region 10mbp
samples="$WORKSPACE_BUCKET/saraj/vntr_samples/passing_sample_ids_v7.1.txt"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN/ACAN_region_50m.bed"
map="$WORKSPACE_BUCKET/saraj/genetic_map/plink.${chr}_w_chr.GRCh38.map"

time python imputation_aou.py \
        --vcf $gt \
        --ref-panel $ref \
        --name sr_ml_filtered_ref_${chr} \
        --window 20 \
	--overlap 2 \
	--chrom $chr \
	--mem 20 \
	--map $map \
	--samples-file $samples \
	--regions-file $regions \
	--batch-size 1000
