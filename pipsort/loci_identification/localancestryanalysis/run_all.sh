#!/bin/bash

for i in {15..20}; do
    echo "Processing iteration $i"
    gsutil cp "${WORKSPACE_BUCKET}/gnomix/outputs/gnomix-chr${i}.msp" ./
    python get_local_ancestry_chrinfo.py $i
    gsutil cp chr${i}* "${WORKSPACE_BUCKET}/pipsort/localancestry/chrinfo/"
    rm gnomix-chr${i}.msp

done
