#!/bin/bash

chr="chr2"
samples="samples/passing_samples_v7.1.csv"
samples_prefix=$(basename $samples | sed 's/.csv//g')
phenotype="cholesterol"
echo "Samples $samples Samples_prefix: $samples_prefix"
summary="summary_gwas_snp_${phenotype}_${samples_prefix}_${chr}.txt"

./aou_gwas.py --phenotype $phenotype \
   --num-pcs 10 \
   --method hail \
   --samples $samples \
   --region chr2 \
   --norm quantile \
   --outdir outputs/${chr}_snp \
   --plot
