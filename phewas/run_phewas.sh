#!/bin/bash

## We use PheTK to run phewas. For more information on installation and documentation see https://github.com/nhgritctran/PheTK/tree/main?tab=readme-ov-file

## To install run:
#pip install PheTK --upgrade

## To run the demo, run:
#python3 -m PheTK.Demo
tr_vcf="../imputation/wdl/data/imputed_chr1.annotated.rh.vcf.gz"
#python run_phewas.py --locus "INS" --chrom "chr11" --tr-vcf $tr_vcf --start 2161569 --n-threads 12

#python run_phewas.py --locus "ACAN" --chrom "chr15" --start 88855424 --n-threads 12 --no-plot --tr-vcf $tr_vcf
for min_case in 10000 5000; do 
time python run_phewas.py --locus "TMCO1" \
                     --chrom "chr1" \
                     --start 165761972 \
                     --n-threads 12 \
                     --no-plot \
                     --min-cases $min_case \
                     --tr-vcf $tr_vcf \
                     --outdir outputs_min_cases_${min_case} \
                     --print-significant-hits
done
#python run_phewas.py --locus "SLC6A3" --chrom "chr5" --start 1393581 --n-threads 12
#python run_phewas.py --locus "EIF3H" --chrom "chr8" --start 116611155 --n-threads 12
#python run_phewas.py --locus "NACA" --chrom "chr12" --start 56717495 --n-threads 12
#python run_phewas.py --locus "CUL4A" --chrom "chr13" --start 113231979 --n-threads 12
#python run_phewas.py --locus "NACA" --chrom "chr12" --start 56717495 --n-threads 12 --min-phecode-count 4
