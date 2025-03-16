#!/bin/bash

time python tr_gwas_aou.py --name sara_tr_gwas \
                           --cohort ALL \
                           --logistic \
                           --pheno-prefix sara_ \
                           --pheno-file ../../../cast-workflows/phewas/phenotypes_min_2000_cases.csv
