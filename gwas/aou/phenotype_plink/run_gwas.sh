#!/bin/bash


# For type 2 diabetes
python convert_phenotype_plink.py \
        --phenotype t1diabetes \
        --is-binary
        #--phenotype t2diabetes \
               
vcf="../../../imputation/wdl/data/sr_ml_filtered_ref_chr11_annotated"
plink2 --pfile $vcf \
	--pheno t1diabetes_pheno_plink.txt \
	--covar t1diabetes_covar_combined.txt \
	--out chr11_t1diabetes_gwas \
	--logistic \
	--covar-variance-standardize
	#--pheno t2diabetes_pheno_plink.txt \
	#--covar t2diabetes_covar_combined.txt \
	#--out chr11_t2diabetes_gwas \

exit 0
# For CAD
python convert_phenotype_plink.py \
        --phenotype CAD
exit 0

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

