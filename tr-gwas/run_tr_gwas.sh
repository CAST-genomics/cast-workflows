#!/bin/bash

#variants="$WORKSPACE_BUCKET/saraj/tr_gwas/variants_DBNDD1.csv"
time python tr_gwas_aou.py --name sara_tr_gwas_eur_male \
                           --cohort EUR \
                           --logistic \
                           --male \
                           --pheno-prefix sara_ \
                           --pheno-file disease_2000_cases_removed_female_specific.csv

                           #--pheno-file disease_phenotypes_significant_hits_2000_cases_all.csv

                           #--phenotype Malignant_neoplasm_of_the_skin
                           #--pheno-file ../../../cast-workflows/phewas/phenotypes_min_2000_cases.csv
                           #--variant-as-covar-file $variants \
                           #--phenotype Atrial_fibrillation,Atrial_fibrillation_and_flutter,Basal_cell_carcinoma,Keratinocyte_carcinoma,Leiomyoma_of_uterus,Schizophrenia,Malignant_neoplasm_of_the_skin
exit 0
time python tr_gwas_aou.py --name sara_tr_gwas \
                           --cohort ALL \
                           --phenotype red_blood_cell_distribution_width,ldl_cholesterol,hdl_cholesterol,bilirubin_direct,mean_corpuscular_hemoglobin,haematocrit,protein,mean_corpuscular_volume,creatinine,platelet_count,alkaline_phosphatase,aspartate_aminotransferase,cholesterol,eosinophil_percent,glucose,mean_platelet_volume,neutrophil_percent,phosphate,triglycerides,urea
