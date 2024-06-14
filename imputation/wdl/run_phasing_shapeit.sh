#!/bin/bash

chr="chr11"
#input=data/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr_AN_AC_re_filled.vcf.gz \
#input="data/aou_${chr}_100samples.vcf.gz"
input="data/100_aou_chr11_imputation.vcf.gz"
# Multiallelic ref
ref="data/${chr}_final_SNP_merged_additional_TRs_AN_AC_filled.vcf.gz"
#ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
# Biallelic ref
#ref="data/${chr}_final_SNP_merged_additional_TRs_biallelic_uniq_id_AN_refilled.vcf.gz"

# To avoid the following error in shapit, run:
#ERROR: AN field is needed in file
#orig_ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
#vcf_base="data/${chr}_final_SNP_merged_additional_TRs_AN_filled.vcf.gz"
#echo "bcftools +fill-tags $orig_ref -Oz -o $vcf_base -- -t AN,AC"
#bcftools +fill-tags $orig_ref -Oz -o $vcf_base -- -t AN
#bgzip -c $vcf_base > $vcf_base.gz
#tabix -p vcf $vcf_base

date
time ./phase_common_static --region ${chr} \
	--input $input \
	--reference $ref \
	--output data/aou_100_${chr}_phased_shapeit_ma.bcf
	#--map genetic_map/genetic_map_hg38_withX.txt.gz
date
