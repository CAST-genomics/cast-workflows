#download targeTR result from google cloud bucket

gsutil cp ${WORKSPACE_BUCKET}/ukb_finemapped_hg38/UKB_finemapped-batch*.filtered.sorted.vcf.gz .
gsutil cp ${WORKSPACE_BUCKET}/eSTR/eSTR_batch1.filtered.sorted.vcf.gz* .


#extract GB field for ukb finemapped str
for i in {1..10}; do { echo -en "CHROM\tPOS\t"; bcftools query -l UKB_finemapped-batch${i}.filtered.sorted.vcf.gz | paste -s;bcftools query -f '%CHROM\t%POS[\t%GB]\n' UKB_finemapped-batch${i}.filtered.sorted.vcf.gz; } > hipstr_batch${i}_GB.txt;done

#extract GB field for eSTR batch1
{ echo -en "CHROM\tPOS\t"; bcftools query -l eSTR_batch1.filtered.sorted.vcf.gz| paste -s;bcftools query -f '%CHROM\t%POS[\t%GB]\n' eSTR_batch1.filtered.sorted.vcf.gz; } > eSTR_GB.txt