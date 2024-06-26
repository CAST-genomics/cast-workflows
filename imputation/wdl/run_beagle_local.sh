#!/bin/bash

# For CBL TR

chr="chr11"
# Have to run this before:
if [ ! -e data/aou_100_${chr}_phased_shapit.vcf.gz ]; then
  bcftools convert -O z data/aou_100_${chr}_phased_shapeit.bcf > data/aou_100_${chr}_phased_shapit.vcf.gz
fi

gt="data/aou_100_${chr}_phased_shapit_ma.vcf.gz"
ref="data/${chr}_final_SNP_merged_additional_TRs.bref3"
out="data/beagle_imputed_shapeit_phased_${chr}_100_aou"

# For ACAN VNTR
#chr="chr15"

##gt="vntr_data/acan_10mbp_aou_185_lrwgs_samples_subset.vcf.gz"
#gt="vntr_data/output_chr15_acan_10mbp_aou_10k_samples.vcf.gz"

##ref="data/${chr}_final_SNP_merged_additional_TRs.bref3"
##ref="vntr_data/phased_ACAN_vntr_snp_650_samples.sorted.vcf.gz"
#ref="vntr_reference/ref_phased_reheader_chr15_acan_vntr.vcf.gz"

#out="data/output_${chr}_10mb_aou_10k_srwgs_w8_o2"

#beagle_5_4_old="beagle.19Apr22.7c0.jar"
beagle_latest="beagle.27May24.118.jar"
beagle=$beagle_latest
date
date > time.txt
java -Xmx8g -jar $beagle \
        gt=$gt \
        ref=$ref \
        out=$out \
        window=8 \
        overlap=2 \
	chrom=$chr
	#map=genetic_map/plink.chr21.GRCh38.map
tabix -p vcf $out.vcf.gz
date >> time.txt
date
