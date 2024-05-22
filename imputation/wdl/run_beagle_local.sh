#!/bin/bash

gt="data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_bialleleic_uniq_id_name.sorted.vcf.gz"
ref="data/chr21_final_SNP_merged_additional_TRs_biallelic_uniq_id_rm_dups.vcf"

java -Xmx12g -jar beagle/beagle.01Mar24.d36.jar \
	gt=$gt \
	ref=$ref \
	out=data/output_chr21_ensemble.vcf \
	window=5 \
	overlap=2
	#map=genetic_map/plink.chr21.GRCh38.map
# 
#gt=data/21_final_SNP_merged_additional_TRs.vcf.gz \
	#map=genetic_map/beagle_chr21_b38.map
