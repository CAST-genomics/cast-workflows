#!/bin/bash

chr="chr11"
new_bref="bref3.01Mar24.d36.jar"
old_bref="bref3.19Apr22.7c0.jar"
bref=$new_bref

zcat data/${chr}_final_SNP_merged_additional_TRs.vcf.gz | \
	java -jar $bref > data/${chr}_final_SNP_merged_additional_TRs.bref3
