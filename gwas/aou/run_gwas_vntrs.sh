#!/bin/bash

echo "" > summary_gwas.txt
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf merged_samples.vcf \
	      --norm zscore \
	      --norm-by-sex \
	      --plot
#./aou_gwas.py --phenotype height \
#	      --num-pcs 10 \
#	      --method hail \
#	      --region chr15:1-99753195\
#	      --norm quantile \
#	      --norm-by-sex \
#	      --plot
#	      #--region chr15:88855424-88857434
#	      #--tr-vcf ${WORKSPACE_BUCKET}/saraj/batch_genotyping/run_8/merged_outputs/run_merge_advntr/d06c3f27-b7b7-41ad-9cfb-8e7161c43fce/call-merge_outputs/merged_samples.vcf

for phenotype in $(tail -n +2 phenotypes_manifest.csv | cut -d, -f1); do
	./aou_gwas.py --phenotype $phenotype \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot
      most_significant_hit=$(tail -n +4 outputs/${phenotype}_associaTR_ALL.gwas.tab | cut -f6 | sort  | awk NF | head -n 1)
      echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> summary_gwas.txt
      echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
done
