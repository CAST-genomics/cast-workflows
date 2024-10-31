#!/bin/bash

python convert_phenotype_plink.py \
        --phenotype t2diabetes \
        --is-binary
               
vcf="../../../imputation/wdl/data/sr_ml_filtered_ref_chr11_annotated"
plink2 --pfile $vcf \
	--pheno t2diabetes_pheno_plink.txt \
	--logistic \
	--covar t2diabetes_covar_combined.txt \
	--out chr11_t2diabetes_gwas \
	--covar-variance-standardize
	#--pheno egfr_pheno_plink.txt \
