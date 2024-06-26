#!/bin/bash

if [ ! -e data/aou_100_${chr}_phased_shapit.vcf.gz ]; then
  bcftools convert -O z data/aou_100_${chr}_phased_shapeit.bcf > data/aou_100_${chr}_phased_shapit.vcf.gz
fi
chr="chr15"
gt="vntr_data/acan_10mbp_aou_185_lrwgs_samples_subset.vcf.gz"
ref="../../../../imputation/cast-workflows/imputation/wdl/vntr_reference/phased_ACAN_vntr_snp_740_samples.sorted.vcf.gz"

out="data/output_${chr}_10mb_aou_185_lrwgs_w8_o2"

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
	chrom=$chr && \
tabix -p vcf $out.vcf.gz
