#!/bin/bash
ref="ref_panel_summary.txt"
echo "chr,snp,tr" | sed 's/,/\t/g' > ${ref}


for (( i=1; i<=22; i++ )); do
    chr="chr$i"
    vcf="${chr}_final_SNP_merged_additional_TRs.vcf.gz"
    echo "counting ${chr}"
    count_snp=$(bcftools view -i 'ID!="."' $vcf |grep -v "^#" | wc -l)
    count_tr=$(bcftools view -i 'ID="."' $vcf |grep -v "^#" | wc -l)

    echo -e "${chr}\t${count_snp}\t${count_tr}" >> ${ref}

done