#!/bin/bash

chr="chr11"
samples="samples/passing_samples_v7.1.csv"
samples_prefix=$(basename $samples | sed 's/.csv//g')
echo "Samples $samples Samples_prefix: $samples_prefix"
summary="summary_gwas_snp_${samples_prefix}_${chr}.txt"

phenotype="type_2_diabetes"
./aou_gwas.py --phenotype $phenotype \
   --num-pcs 10 \
   --method hail \
   --samples $samples \
   --norm quantile \
   --outdir outputs/${chr}_snp \
   --plot
