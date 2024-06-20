#!/bin/bash

# Have to run this before:
#bcftools convert -O z data/aou_100_chr11_phased_shapeit_ma_.bcf > data/aou_100_chr11_phased_shapit_ma.vcf.gz
chr="chr15"
gt="vntr_data/acan_10mbp_aou_185_lrwgs_samples_subset.vcf.gz"
ref="../../../../imputation/cast-workflows/imputation/wdl/vntr_reference/phased_ACAN_vntr_snp_740_samples.sorted.vcf.gz"

out="data/output_${chr}_10mb_aou_185_lrwgs_w8_o2"

beagle_5_4_old="beagle.19Apr22.7c0.jar"
beagle=$beagle_5_4_old
java -Xmx4g -jar $beagle \
        gt=$gt \
        ref=$ref \
        out=$out \
        window=8 \
        overlap=2 \
	chrom=$chr && \
tabix -p vcf $out.vcf.gz
