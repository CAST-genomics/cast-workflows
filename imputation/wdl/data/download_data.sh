#!/bin/bash

for chr in "chr15"; do

  gt_file="ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  gt_file_chr="ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr.vcf"
  gt_url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/${gt_file}"
  
  if [ ! -e $gt_file ]; then
    wget $gt_url
    wget $gt_url.tbi
  fi
  
  if [ ! -e $gt_file_chr ]; then
   zcat $gt_file | grep "^#" | sed 's/contig=<ID=15>/contig=<ID=chr15>/g' > $gt_file_chr
   zcat $gt_file | grep -v "^#" | awk '{print("chr"$0)}' >> $gt_file_chr
  fi
  
  bgzip -c $gt_file_chr > $gt_file_chr.gz
  tabix -p vcf $gt_file_chr.gz
  rm $gt_file_chr
  
exit 0

  ref_file="${chr}_final_SNP_merged_additional_TRs.vcf.gz"
  ref_url="https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/$ref_file"
  
  if [ ! -e $ref_file ]; then
    wget $ref_url
    wget $ref_url.tbi
  fi
done
