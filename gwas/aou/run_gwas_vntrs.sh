#!/bin/bash



./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf merged_samples.vcf \
	      --norm zscore \
	      --norm-by-sex \
	      --plot
	      #--tr-vcf ${WORKSPACE_BUCKET}/saraj/batch_genotyping/run_8/merged_outputs/run_merge_advntr/d06c3f27-b7b7-41ad-9cfb-8e7161c43fce/call-merge_outputs/merged_samples.vcf

./aou_gwas.py --phenotype glucose \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot

./aou_gwas.py --phenotype calcium \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot

./aou_gwas.py --phenotype ALT \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot


./aou_gwas.py --phenotype platelet_count \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot
	      

./aou_gwas.py --phenotype ldl_cholesterol \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot
	      
./aou_gwas.py --phenotype hdl_cholesterol \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot

./aou_gwas.py --phenotype cholesterol \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot
exit 0

./aou_gwas.py --phenotype colon_polyp \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf merged_samples.vcf \
	      --norm quantile \
	      --plot
