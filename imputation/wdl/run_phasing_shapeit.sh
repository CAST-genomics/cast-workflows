#!/bin/bash

#input=data/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr_AN_AC_re_filled.vcf.gz \
input="data/aou_chr21_100samples.vcf.gz"
# Multiallelic ref
ref="data/chr21_final_SNP_merged_additional_TRs_AN_AC_re_filled.vcf.gz"
# Biallelic ref
#ref="data/chr21_final_SNP_merged_additional_TRs_biallelic_uniq_id_AN_refilled.vcf.gz"

# To avoid the following error in shapit, run:
#ERROR: AN field is needed in file
#biallelic_ref="data/chr21_final_SNP_merged_additional_TRs_biallelic_uniq_id.vcf.gz"
#vcf_base="data/chr21_final_SNP_merged_additional_TRs_biallelic_uniq_id_AN_refilled.vcf"
#echo "bcftools +fill-tags $biallelic_ref -Oz -o $vcf_base -- -t AN,AC"
#bcftools +fill-tags $biallelic_ref -Oz -o $vcf_base -- -t AN,AC
#bgzip -c $vcf_base > $vcf_base.gz
#tabix -p vcf $vcf_base.gz

date
time ./phase_common_static --region chr21 \
	--input $input \
	--reference $ref \
	--output data/aou_100_chr21_phased_shapit_ma.bcf
	#--map genetic_map/genetic_map_hg38_withX.txt.gz
date
