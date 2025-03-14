#!/bin/bash

while read line; do
    phecode=$(echo $line | awk -F, '{print($1)}')
    string=$(echo $line | awk -F, '{print($4)}')
    string="sara_$string"
    if [ ! -e phenocovars/${string}_phenocovar.csv ]; then
        echo "phecode $phecode string $string"
    fi
    ./aou_phenotype_preprocessing.py \
            --phenotype $string \
            --concept-id $phecode \
            --samples passing_samples_v7.1.csv \
            --phecode-cohort-file ../../../phewas/my_phecode_counts.csv \
            --cohort-file ../../../phewas/cohort_with_covariates.csv \
            --phecode-count 2 \
            --outdir phenocovars \
            --verbose
done < <(tail -n +2 ../../../phewas/phenotypes_min_1000_cases.csv)
