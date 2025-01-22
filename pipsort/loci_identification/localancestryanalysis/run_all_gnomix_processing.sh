#!/bin/bash

for i in {2..2}; do
    echo "Processing iteration $i"
    gsutil cp "${WORKSPACE_BUCKET}/gnomix/outputs/gnomix-chr${i}.msp" ./
    python get_local_ancestry_chrinfo_per_cohort.py $i NOT_AFR_BLACK.csv
    python get_local_ancestry_chrinfo_per_cohort.py $i AFR_BLACK.csv
    gsutil cp chr${i}*AFR_BLACK.csv "${WORKSPACE_BUCKET}/pipsort/localancestry/chrinfo/"
    rm gnomix-chr${i}.msp

done
