#!/bin/bash

chr="chr21"
panel_base="data/${chr}_final_SNP_merged_additional_TRs"
panel_vcf="$panel_base.vcf.gz"
panel_vcf_uniq_id="${panel_base}_biallelic_uniq_id.vcf.gz"
panel_vcf_biallelic="${panel_vcf_base}_biallelic.vcf"
panel_vcf_uniq_id_sorted="${panel_vcf_base}_bialleleic_uniq_id.sorted.vcf.gz"

gt_vcf_base="data/ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr"
gt_vcf="$gt_vcf_base.vcf.gz"
gt_vcf_biallelic="${gt_vcf_base}_biallelic.vcf"
gt_vcf_uniq_id="${gt_vcf_base}_bialleleic_uniq_id.vcf.gz"
gt_vcf_uniq_id_sorted="${gt_vcf_base}_bialleleic_uniq_id.sorted.vcf.gz"

vcf=$panel_vcf
vcf_biallelic=$panel_vcf_biallelic
vcf_uniq_id=$panel_vcf_uniq_id
vcf_uniq_id_sorted=$panel_vcf_uniq_id_sorted

echo "Changing multiallelic loci to biallelic"
bcftools norm -m-any \
		$vcf \
		> $vcf_biallelic

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
