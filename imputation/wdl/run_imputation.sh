#!/bin/bash

chr="chr15"
# EnsembleTR ref
#ref="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/${chr}_final_SNP_merged_additional_TRs.vcf.gz"

# PVNTR ref
#ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ref_phased_output_chr15_acan_vntr_apr_22.bref3"
ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz"
#ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/phased_ACAN_vntr_snp_650_samples.sorted.vcf.gz"

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
gt="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.${chr}.vcf.bgz"
#gt="gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/saraj/acan_data_all/output_chr15_acan_50mbp_aou_10_samples.vcf.gz"
# VCF file with 100 samples
#gt="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/aou_data_subset/100_aou_${chr}_imputation.vcf.gz"
# VCF file for all aou samples, only acan region 10mbp
#gt="$WORKSPACE_BUCKET/saraj/vntr_samples/aou_all_chr15_acan_10mbp.vcf.gz"
#samples="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/subset_samples/aou_subset_samples_100.txt"
samples="$WORKSPACE_BUCKET/saraj/vntr_samples/sample_ids.txt"
#samples="$WORKSPACE_BUCKET/saraj/vntr_samples/vntr_samples_275_test.txt"
#regions="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/CBL_extended_region.bed"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_50m.bed"

time python imputation_aou.py \
        --vcf $gt \
        --ref-panel $ref \
        --name output_${chr}_acan_50mbp_aou_all_srwgs_samples \
        --window 10 \
	--overlap 2 \
	--chrom $chr \
	--mem 200 \
	--samples-file $samples \
	--regions-file $regions
#	--cromwell

