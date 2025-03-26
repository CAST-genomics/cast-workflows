#!/bin/bash

#TODO: Add Malignant_neoplasm_of_the_pancreas

time python tr_gwas_aou.py --name sara_tr_gwas \
                           --cohort EUR,AFR \
                           --logistic \
                           --pheno-prefix sara_ \
                           --phenotype Atrial_fibrillation,Atrial_fibrillation_and_flutter,Basal_cell_carcinoma,Keratinocyte_carcinoma,Leiomyoma_of_uterus,Schizophrenia,Malignant_neoplasm_of_the_skin
                           #--pheno-file ../../../cast-workflows/phewas/phenotypes_min_2000_cases.csv
