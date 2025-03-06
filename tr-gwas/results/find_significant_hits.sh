#!/bin/bash

# Get lines with significant association.
grep ADD *.tab | sort -k16 -g | grep -v "\sNA\s" | awk -v thr=1e-6 '$16 < thr' > significant_hits.tab

# If p_g_vntr_ids.txt is not present, run
# bcftools view ../../../../cast-workflows/imputation/wdl/data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz | grep -v "^##"| awk '{print $3}' > p_g_vntr_ids.txt

# Find the VNTRs among the significant hits
cat significant_hits.tab | awk '{print $3}' | grep -f p_g_vntr_ids.txt > vids_with_significant_hits.txt
grep -f vids_with_significant_hits.txt significant_hits.tab > significant_hits_vntrs.tab
