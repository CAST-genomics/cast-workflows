#!/bin/bash

# Get lines with significant association.
grep ADD *.tab | sort -k15 -g | grep -v "\sNA\s" | awk -v thr=5e-8 '$15 < thr' | \
        sed 's/_gwas\.tab:[0-9]*//g' | \
        sed 's/sara_//g' | \
        sed 's/_EUR_WHITE//g' | \
        sed 's/_AFR_BLACK//g' | \
        sed 's/_AMR_HISPANIC//g' | \
        sed 's/_passing_samples_v7.1//g' > significant_hits.tab > significant_hits.tab

echo "Number of unique phenotypes"
cat significant_hits.tab | awk '{print $1}' | awk -F: '{print $1}'| sort | uniq -c | sort -rg | wc -l
cat significant_hits.tab | awk '{print $1}' | awk -F: '{print $1}'| sort | uniq -c | sort -rg  > phenotypes_and_counts_in_significant_hits.txt

echo "Number of unique loci"
cat significant_hits.tab | awk '{print $3}' | sort | uniq -c | sort -rg | wc -l
cat significant_hits.tab | awk '{print $3}' | sort | uniq -c | sort -rg > loci_and_counts_in_significant_hits.txt
