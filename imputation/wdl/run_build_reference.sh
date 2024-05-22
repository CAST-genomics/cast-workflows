#!/bin/bash

vcf="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/21_final_SNP_merged_additional_TRs.vcf.gz"
ref="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
map="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/genetic_map_hg38_withX.txt.gz"
python build_reference_aou.py \
	--name tr_build_reference \
	--vcf $vcf \
	--ref-panel $ref \
	--map $map
