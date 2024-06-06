#!/bin/bash



beagle_vcf="data/output_chr11_aou_100_phasing_impute_beagle_apr_22_bref.vcf.gz"
beagle_cbl="data/CBL_aou_100_phasing_impute_beagle_apr_22_bref.vcf.gz"
shapeit_vcf="data/output_chr11_aou_100_phasing_impute_beagle_apr_22_bref.vcf.gz"
shapeit_cbl="data/CBL_aou_100_phase_shapeit_impute_beagle_apr_22_bref.vcf.gz"
#bcftools view -O z -r chr11:119206290 $beagle_vcf > $beagle_cbl
tabix -p vcf $beagle_cbl
#bcftools view -O z -r chr11:119206290 $shapeit_vcf > $shapeit_cbl
tabix -p vcf $shapeit_cbl

comparestr="../../../../trtools/TRTools/trtools/compareSTR/compareSTR.py"
vcf1="data/CBL_test.filtered.sorted.vcf.gz"
vcf2=$beagle_cbl
out="data/compare_gt/CBL_compare_beagle_hipstr"


python $comparestr \
	--vcf1 $vcf1 --vcf2 $vcf2 --out $out
