#!/bin/bash

#variants="$WORKSPACE_BUCKET/saraj/tr_gwas/variants_DBNDD1.csv"
blood_phenotypes="red_blood_cell_distribution_width,ldl_cholesterol,hdl_cholesterol,bilirubin_direct,mean_corpuscular_hemoglobin,haematocrit,protein,mean_corpuscular_volume,creatinine,platelet_count,alkaline_phosphatase,aspartate_aminotransferase,cholesterol,eosinophil_percent,glucose,mean_platelet_volume,neutrophil_percent,phosphate,triglycerides,urea"
pheno_file="disease_phenotypes_significant_hits_2000_cases_all.csv"
time python tr_gwas_aou.py --name sara_tr_gwas_amr_male \
                           --cohort AMR \
                           --logistic \
                           --male \
                           --pheno-prefix sara_ \
                           --pheno-file disease_2000_cases_removed_female_specific_error_removed.csv
                           #--variant-as-covar-file $variants \
