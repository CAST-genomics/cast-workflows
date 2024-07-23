#!/bin/bash

chr="chr11"
ref="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/${chr}_final_SNP_merged_additional_TRs.vcf.gz"

## Run with 1000 genomes data
#gt="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/query/ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr.vcf.gz"
#samples="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/1000genome_10sample.txt"
#python imputation_aou.py \
#	--name tr_imputation \
#	--vcf $vcf \
#	--ref-panel $ref \
#	--samples-file $samples \
#	--regions-file $regions

## Run with AoU data
# VCF file with all samples
#gt="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.${chr}.vcf.bgz"
# VCF file with 100 samples
gt="gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/tr_imputation/tr_imputation/aou_data_subset/100_aou_${chr}_imputation.vcf.gz"
samples="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/aou_imputation_5sample.txt"
regions="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/CBL_extended_region.bed"

python imputation_aou.py \
        --vcf $gt \
        --ref-panel $ref \
        --name output_${chr}_ensemble \
        --window 5 \
	--samples-file $samples \
	--regions-file $regions

