#!/bin/bash

#aou_220k_imputed="data/merged_samples_all_rh.sorted.vcf.gz"
#aou_220k_imputed="../../imputation/wdl/data/imputed_p_g_vntrs_srwgs_sr_ml_filter.sorted.annotated.vcf.gz"
aou_220k_imputed="../../imputation/wdl/data/chr_11_imputed.sorted.annotate.vcf.gz"

#aou_220k_imputed="data/vntr_ref.sorted.vcf.gz"
#ref="data/ACAN_ref_lrwgs.sorted.vcf.gz"
ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"


samples="samples/passing_samples_v7.1.csv"
#samples="samples/EUR_WHITE.csv"
#samples="samples/AFR_BLACK.csv"

samples_prefix=$(basename $samples | sed 's/.csv//g')
echo "Samples $samples Samples_prefix: $samples_prefix"

echo "" > summary_gwas.txt
for phenotype in $(tail -n +2 phenotypes_manifest.csv  | cut -d, -f1); do
     echo "Running gwas for $phenotype"
     ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --norm quantile \
           --is-imputed \
           --plot
     most_significant_hit=$(tail -n +4 outputs/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)
     echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> summary_gwas.txt
     echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
done

exit 0

# Running associatr
echo "running gwas for imputed calls"
samples="samples/AFR_BLACK.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --samples $samples \
	      --is-imputed \
	      --norm-by-sex \
	      --norm quantile \
	      --annotations annotations_acan.txt \
	      --plot
exit 0
samples="samples/EUR_WHITE.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --annotations annotations_acan.txt \
	      --samples $samples \
	      --is-imputed \
	      --norm-by-sex \
	      --norm quantile \
	      --plot
samples="samples/passing_samples_v7.1.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --annotations annotations_acan.txt \
	      --samples $samples \
	      --is-imputed \
	      --norm-by-sex \
	      --norm quantile \
	      --plot
echo "running gwas for advntr calls"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $lrwgs_data \
	      --annotations annotations_acan.txt \
	      --norm quantile \
	      --samples $samples \
	      --norm-by-sex \
	      --annotations annotations_acan.txt \
	      --plot
#	      --tr-vcf merged_samples.vcf \
exit 0
echo "running gwas for reference"
samples="samples/passing_samples_v7.1.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $ref\
	      --samples $samples \
	      --norm-by-sex \
	      --norm quantile \
	      --annotations annotations_acan.txt \
	      --plot
exit 0
