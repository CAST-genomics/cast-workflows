#!/bin/bash

#chr="chr11"
chr="chr14"
#input=data/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr_AN_AC_re_filled.vcf.gz \
#input="data/aou_${chr}_100samples.vcf.gz"
#input="data/100_aou_${chr}_imputation.vcf.gz"
# Hipstr calls
#input="data/UKB_finemapped-batch3.filtered.sorted.vcf.gz"
input="data/output_chr14_tr_10mbp_aou_100.vcf.gz"
# Multiallelic ref
ref="data/${chr}_final_SNP_merged_additional_TRs_AN_AC_filled.vcf.gz"
#ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
# Biallelic ref
#ref="data/${chr}_final_SNP_merged_additional_TRs_biallelic_uniq_id_AN_refilled.vcf.gz"

# To avoid the following error in shapit, run:
#ERROR: AN field is needed in file
#orig_ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
#vcf_base="data/${chr}_final_SNP_merged_additional_TRs_AN_AC_filled.vcf.gz"
#echo "bcftools +fill-tags $orig_ref -Oz -o $vcf_base -- -t AN,AC"
#bcftools +fill-tags $orig_ref -Oz -o $vcf_base -- -t AN,AC
##bgzip -c $vcf_base > $vcf_base.gz
#tabix -p vcf $vcf_base

output="data/aou_all_${chr}_phased_shapeit.bcf"
date
time ./phase_common_static --region ${chr} \
	--input $input \
	--reference $ref \
	--output $output
	#--map genetic_map/genetic_map_hg38_withX.txt.gz
date
