#!/bin/bash

chr="chr13"
#aou_220k_imputed="../../imputation/wdl/data/imputed_chr15_p_g_vntrs_srwgs_sr_ml_filter.sorted.vcf.gz"
#aou_220k_imputed="../../imputation/wdl/data/imputed_chr11.annotated.rh.vcf.gz"
#aou_220k_imputed="../../imputation/wdl/data/imputed_samples_chr13.sorted.annotate.rh.vcf.gz"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"

ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"

#for samples in "samples/AFR_BLACK.csv" "samples/EUR_WHITE.csv"; do
  samples="samples/passing_samples_v7.1.csv"
  #samples="samples/AFR_BLACK.csv"
  #samples="samples/EUR_WHITE.csv"
  samples_prefix=$(basename $samples | sed 's/.csv//g')
  echo "Samples $samples Samples_prefix: $samples_prefix"
  summary="summary_gwas_${samples_prefix}_${chr}_v2.txt"


  echo "" >> $summary
  #for phenotype in $(tail -n +2 phenotypes_manifest.csv  | cut -d, -f1); do
     phenotype="red_blood_cell_distribution_width"
     #phenotype="glaucoma"
     #phenotype="type_2_diabetes"
     #phenotype="height"
     #snp_gwas_file="data/all_by_all/df_dump_${chr}_EM_202.2.csv"
     snp_gwas_file="data/all_by_all/df_dump_${chr}_${phenotype}.csv"
    
     echo "Running gwas for $phenotype"
     ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --norm quantile \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --snp-gwas-file $snp_gwas_file \
           --plot
           #--annotations annotations_acan.txt \
     most_significant_hit=$(tail -n +4 outputs/${chr}/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)
     echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> $summary
     echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
#done
