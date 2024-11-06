#!/bin/bash

vcf="data/imputed_chr15_p_g_vntrs_srwgs_sr_ml_filter.sorted.vcf.gz"
ref_vcf="data/vntr_ref_chr15.sorted.annotated.vcf.gz"
out_prefix="imputed_chr15"

python ../../../TRTools/trtools/annotaTR/annotaTR.py --vcf ${vcf} \
                 --ref-panel ${ref_vcf} \
                 --out ${out_prefix}_annotated \
                 --vcftype advntr \
                 --outtype pgen \
                 --dosages bestguess_norm

