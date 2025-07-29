#!/bin/bash

phen=$1
loci=$2

while IFS= read -r line; do
    plink2 \
        --bfile "$line" \
        --make-pgen \
        --out "$line"
done < "$loci"
