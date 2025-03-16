#!/bin/bash

# Get lines with significant association.
grep ADD *.tab | sort -k16 -g | grep -v "\sNA\s" | awk -v thr=5e-8 '$16 < thr' > significant_hits.tab

echo "Number of unique phenotypes"
cat significant_hits.tab | awk '{print $1}' | awk -F: '{print $1}'| sort | uniq -c | sort -rg | wc -l

echo "Number of unique loci"
cat significant_hits.tab | awk '{print $3}' | sort | uniq -c | sort -rg | wc -l
