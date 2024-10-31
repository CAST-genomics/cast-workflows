#!/bin/bash

# For CAD

python convert_phenotype_plink.py \
        --phenotype CAD


# For height
python convert_phenotype_plink.py \
        --phenotype height
               
vcf="../../../imputation/wdl/data/sr_ml_filtered_ref_chr15_annotated"
plink2 --pfile $vcf \
	--pheno height_pheno_plink.txt \
	--linear \
	--covar height_covar_combined.txt \
	--out chr15_height_gwas \
	--covar-variance-standardize

exit 0
# For type 2 diabetes
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
