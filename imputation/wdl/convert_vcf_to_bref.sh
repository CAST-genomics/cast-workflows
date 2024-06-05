#!/bin/bash

chr="chr15"
new_bref="bref3.01Mar24.d36.jar"
old_bref="bref3.19Apr22.7c0.jar"
bref=$new_bref


#vcf_base="data/${chr}_final_SNP_merged_additional_TRs"
#vcf_base="vntr_reference/ACAN_vntr_snp.sorted"
vcf_base="vntr_reference/ref_phased_output_chr15_acan_vntr_apr_22"

zcat $vcf_base.vcf.gz | \
     java -jar $bref > $vcf_base.bref3
