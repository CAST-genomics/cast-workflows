#!/bin/bash

chr="chr1"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"

ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"

for samples in "samples/passing_samples_v7.1.csv"; do
  #samples="samples/passing_samples_v7.1.csv"
  #samples="samples/AFR_BLACK.csv"
  #samples="samples/EUR_WHITE.csv"
  samples_prefix=$(basename $samples | sed 's/.csv//g')
  echo "Samples $samples Samples_prefix: $samples_prefix"
  summary="summary_gwas_${samples_prefix}_${chr}_v2.txt"


  echo "" >> $summary
     phenotype="glaucoma"
     snp_gwas_file="data/all_by_all/df_dump_${chr}_${phenotype}.csv"
    
     echo "Running gwas for $phenotype"
     ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --plot-genotype-phenotype \
           --downsample \
           --merge-bins \
           --violin \
           --plot
      #     --snp-gwas-file $snp_gwas_file \
           #--norm quantile \
     most_significant_hit=$(tail -n +4 outputs/${chr}/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)
     echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> $summary
     echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
done
