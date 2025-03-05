#!/bin/bash

python aou_phenotype_preprocessing.py \
	--phenotype glaucoma \
    --concept-id SO_375.1 \
    --phecode-cohort-file ../../../phewas/my_phecode_counts.csv \
    --phecode-count 2 \
    --verbose
