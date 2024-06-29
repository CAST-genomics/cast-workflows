#!/bin/bash

# For CBL TR

chr="chr14"
# Have to run this before:
gt_base="data/aou_100_${chr}_phased_shapeit"
if [ ! -e "$gt_base.vcf.gz" ]; then
  bcftools convert -Oz "$gt_base.bcf" > "$gt_base.vcf"
  tabix -p vcf "$gt_base.vcf.gz"
fi

gt="data/aou_100_${chr}_phased_shapeit.vcf.gz"
#gt="../../../../imputation/cast-workflows/imputation/wdl/data/aou_100_${chr}_phased_shapit_ma.vcf.gz"
#ref="../../../../imputation/cast-workflows/imputation/wdl/data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
#ref="data/${chr}_final_SNP_merged_additional_TRs.bref3"
#ref="data/${chr}_final_SNP_merged_additional_TRs_AN_AC_filled.vcf.gz"
ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
out="data/beagle_imputed_shapeit_phased_${chr}_100_aou_test"

gt="data/output_chr14_tr_50mbp_aou_100.vcf.gz"
out="data/beagle_imputed_and_phased_${chr}_100_aou_50mbp"


# For ACAN VNTR
#chr="chr15"

##gt="vntr_data/acan_10mbp_aou_185_lrwgs_samples_subset.vcf.gz"
#gt="vntr_data/output_chr15_acan_10mbp_aou_10k_samples.vcf.gz"

##ref="data/${chr}_final_SNP_merged_additional_TRs.bref3"
##ref="vntr_data/phased_ACAN_vntr_snp_650_samples.sorted.vcf.gz"
#ref="vntr_reference/ref_phased_reheader_chr15_acan_vntr.vcf.gz"

#out="data/output_${chr}_10mb_aou_10k_srwgs_w8_o2"

beagle_5_4_old="beagle.19Apr22.7c0.jar"
beagle_latest="beagle.27May24.118.jar"
beagle=$beagle_latest
date
date > time.txt
java -Xmx8g -jar $beagle \
        gt=$gt \
        ref=$ref \
        out=$out \
        window=5 \
        overlap=2 \
	chrom=$chr
	#map=genetic_map/plink.chr21.GRCh38.map
tabix -p vcf $out.vcf.gz
date >> time.txt
date
