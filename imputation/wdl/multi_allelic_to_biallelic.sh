#!/bin/bash

panel_vcf="data/chr21_final_SNP_merged_additional_TRs_biallelic.vcf"
panel_vcf_uniq_id="data/chr21_final_SNP_merged_additional_TRs_biallelic_uniq_id.vcf.gz"
gt_vcf="data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
gt_vcf_biallelic="data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_biallelic.vcf"
gt_vcf_uniq_id="data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_bialleleic_uniq_id.vcf.gz"
gt_vcf_uniq_id_sorted="data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_bialleleic_uniq_id.sorted.vcf.gz"
vcf=$gt_vcf
vcf_biallelic=$gt_vcf_biallelic
vcf_uniq_id=$gt_vcf_uniq_id
vcf_uniq_id_sorted=$gt_vcf_uniq_id_sorted

echo "Changing multiallelic loci to biallelic"
bcftools norm -m-any \
		$gt_vcf \
		> $gt_vcf_biallelic

echo "Assigning uniq ids"

for CHR in {21}; do
    bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    $vcf_biallelic  \
    -Oz -o $vcf_uniq_id
done

echo "sort and index"
zcat $vcf_uniq_id | \
	vcf-sort | \
	bgzip -c > \
	$vcf_uniq_id_sorted && \
	tabix -p vcf $vcf_uniq_id_sorted
