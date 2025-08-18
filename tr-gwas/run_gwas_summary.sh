#!/bin/bash

phecode_file="../../../cast-workflows/phewas/phenotypes_min_2000_cases.csv"
#summary="results_v8_eur_female_significant_in_v2/significant_hits_removed_male_specific_phenos.tab"
#summary="results_v9_eur_male_significant_in_v2/significant_hits_filtered_removed_female_specific.tab"
cohort="AMR"
summary="results_v7_amr_significant_in_v2/significant_hits_filtered.tab"
python gwas_summary_plot.py --gwas-summary $summary \
                            --cohort $cohort \
                            --phecode-to-name-file $phecode_file

cohort="AFR"
summary="results_v6_afr_significant_in_v2/significant_hits_filtered.tab"
python gwas_summary_plot.py --gwas-summary $summary \
                            --cohort $cohort \
                            --phecode-to-name-file $phecode_file

cohort="EUR"
summary="results_v3_eur_significant_in_v2/significant_hits_filtered.tab"
python gwas_summary_plot.py --gwas-summary $summary \
                            --cohort $cohort \
                            --phecode-to-name-file $phecode_file


cohort="ALL"
summary="results_v4_blood/significant_hits_filtered.tab"
python gwas_summary_plot.py --gwas-summary $summary \
                            --cohort $cohort \
                            --is-continuous \
                            --phecode-to-name-file $phecode_file

