#!/bin/bash

chr="chr16"
samples="samples/EUR_WHITE.csv"
samples_prefix=$(basename $samples | sed 's/.csv//g')
phenotype="Malignant_neoplasm_of_the_skin"
echo "Samples $samples Samples_prefix: $samples_prefix"
summary="summary_gwas_snp_${phenotype}_${samples_prefix}_${chr}.txt"
tr_vcf="data/imputed_chr16.sorted.annotated.vcf.gz"
./aou_gwas.py --phenotype $phenotype \
   --num-pcs 10 \
   --method hail \
   --tr-vcf $tr_vcf \
   --vntr-as-covariate chr16_90007850_DBNDD1 \
   --samples $samples \
   --region chr16:89994871-90029456 \
   --norm quantile \
   --outdir outputs/${chr}_snp \
   --plot
