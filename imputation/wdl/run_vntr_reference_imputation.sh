#!/bin/bash


chr="chr15"
#gt="data/aou_${chr}_100samples.vcf.gz"
gt="data/ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr.vcf.gz"

#ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
#ref="vntr_reference/ACAN_vntr_snp.sorted.vcf.gz"
ref="vntr_reference/ref_phased_output_chr15_acan_vntr_apr_22.bref3"
out="data/output_${chr}_acan_vntr_apr_22"

beagle_5_4="beagle.01Mar24.d36.jar"
beagle_5_4_old="beagle.19Apr22.7c0.jar"
beagle_5_3="beagle.08Feb22.fa4.jar"
beagle=$beagle_5_4_old
date
date > time.txt
time java -Xmx10g -jar $beagle \
	ref=$ref \
	gt=$gt \
	out=$out \
	window=10 \
	chrom=chr15 \
	overlap=2
	#ref=$ref \
tabix -p vcf $out.vcf.gz
