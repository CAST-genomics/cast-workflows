#!/bin/bash

mkdir samples genotypes data phenocovars outputs
gsutil cp $WORKSPACE_BUCKET/samples/EUR_WHITE.csv samples/
gsutil cp $WORKSPACE_BUCKET/samples/passing_samples_v7.1.csv samples/
gsutil cp $WORKSPACE_BUCKET/phenotypes/case/sara_Malignant_neoplasm_of_the_skin_phenocovar.csv phenocovars/
mv phenocovars/sara_Malignant_neoplasm_of_the_skin_phenocovar.csv phenocovars/Malignant_neoplasm_of_the_skin_phenocovar.csv
gsutil cp $WORKSPACE_BUCKET/saraj/imputed_genotypes_dump/* genotypes/
gsutil cp $WORKSPACE_BUCKET/saraj/imputation_output/chr16/imputed_chr16.sorted.annotated.vcf.gz* data/
