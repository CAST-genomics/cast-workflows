#!/bin/bash

chr="chr15"
ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/p_g_vntrs/phased/chr15/vntr_ref.sorted.vcf.gz"


# VCF file with all samples
gt="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.${chr}.vcf.bgz"
samples="$WORKSPACE_BUCKET/saraj/vntr_samples/passing_sample_ids_v7.1.txt"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN/ACAN_region_50m.bed"
map="$WORKSPACE_BUCKET/saraj/genetic_map/plink.chr15_w_chr.GRCh38.map"

time python imputation_aou.py \
        --vcf $gt \
        --ref-panel $ref \
        --name sr_ml_filtered_ref \
        --window 20 \
	--overlap 2 \
	--chrom $chr \
	--vid 290964 \
	--motif CTGCCCCTGGAGTAGAGGACATCAGCGGGCTTCCTTCTGGAGAAGTTCTAGAGACTG \
	--mem 20 \
	--map $map \
	--samples-file $samples \
	--regions-file $regions \
	--batch-size 1000


