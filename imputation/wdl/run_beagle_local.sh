#!/bin/bash

chr="chr21"
#gt="data/ALL.${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_bialleleic_uniq_id_name.sorted.vcf.gz"
gt="data/ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr.vcf.gz"
ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"

beagle_5_4="beagle.01Mar24.d36.jar"
beagle_5_4_old="beagle.19Apr22.7c0.jar"
beagle_5_3="beagle.08Feb22.fa4.jar"
beagle=$beagle_5_4_old
java -Xmx10g -jar $beagle \
	gt=$gt \
	ref=$ref \
	out=data/output_${chr}_ensemble \
	window=5 \
	overlap=2
	#overlap=2
	#map=genetic_map/plink.chr21.GRCh38.map
# 
#gt=data/21_final_SNP_merged_additional_TRs.vcf.gz \
	#map=genetic_map/beagle_chr21_b38.map
