#!/bin/bash

vcf="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_bialleleic_uniq_id.sorted.vcf.gz"
ref="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/chr21_final_SNP_merged_additional_TRs_biallelic_uniq_id.sorted.vcf.gz"
python imputation_aou.py \
	--name tr_imputation \
	--vcf $vcf \
	--ref-panel $ref
