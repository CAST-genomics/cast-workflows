#!/bin/bash

#aou_220k_imputed="../../imputation/wdl/data/imputed_p_g_vntrs_srwgs_sr_ml_filter.sorted.annotated.vcf.gz"
#aou_220k_imputed="../../imputation/wdl/data/chr_11_imputed_rh.sorted.annotate.vcf.gz"
aou_220k_imputed="../../imputation/wdl/data/imputed_samples_chr13.sorted.annotate.rh.vcf.gz"
chr="chr13"

ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"

#for samples in "samples/passing_samples_v7.1.csv" "samples/AFR_BLACK.csv"; do
  samples="samples/passing_samples_v7.1.csv"
  samples_prefix=$(basename $samples | sed 's/.csv//g')
  echo "Samples $samples Samples_prefix: $samples_prefix"
  summary="summary_gwas_${samples_prefix}_${chr}.txt"


  echo "" > $summary
  for phenotype in $(tail -n +2 phenotypes_manifest.csv  | cut -d, -f1); do
     #phenotype="red_blood_cell_distribution_width"
     #phenotype="hemoglobin_a1c"
     #phenotype="diabetes"
     #phenotype="glucose"
     #snp_gwas_file="data/all_by_all/df_dump_${chr}_3019897.csv"
     
     #phenotype="haematocrit"
     #snp_gwas_file="data/all_by_all/df_dump_${chr}_3023314.csv"

     #phenotype="mean_platelet_volume"
     #snp_gwas_file="data/all_by_all/df_dump_${chr}_3043111.csv"
     
     #phenotype="alkaline_phosphatase"
     #snp_gwas_file="data/all_by_all/df_dump_${chr}_3035995.csv"
     
     #phenotype="urea"
     #snp_gwas_file="data/all_by_all/df_dump_${chr}_3013682.csv"

     echo "Running gwas for $phenotype"
     ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --norm quantile \
           --is-imputed \
           --outdir outputs/${chr} \
           --plot
     #      --snp-gwas-file $snp_gwas_file \
     most_significant_hit=$(tail -n +4 outputs/${chr}/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)
     echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> $summary
     echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
  done
#done
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
