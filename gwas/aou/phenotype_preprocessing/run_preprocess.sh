#!/bin/bash

while read line; do
    phecode=$(echo $line | awk -F, '{print($1)}')
    # Replace spaces in the name with underlines.
    # Remove square brackets in the name if present.
    string=$(echo $line | awk -F, '{print($4)}' | sed 's/ /_/g' | sed 's/\[//g' | sed 's/\]//g')
    echo "phecode $phecode string $string"
    ./aou_phenotype_preprocessing.py \
            --phenotype $string \
            --concept-id $phecode \
            --phecode-cohort-file ../../../phewas/my_phecode_counts.csv \
            --phecode-count 2 \
            --outdir phenocovars \
            --verbose
done < <(tail -n +2 ../../../phewas/common_phenotypes_10k_cases.csv)
