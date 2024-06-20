#!/bin/bash

chr="chr15"
ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz"
#ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/phased_ACAN_vntr_snp_650_samples.sorted.vcf.gz"

# VCF file with all samples
gt="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.${chr}.vcf.bgz"
samples="$WORKSPACE_BUCKET/saraj/vntr_samples/sample_ids_aou_10000.txt"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_10m.bed"

time python imputation_aou.py \
        --vcf $gt \
        --ref-panel $ref \
        --name output_${chr}_acan_10mbp_aou_10k_samples \
        --window 5 \
	--overlap 2 \
	--chrom chr15 \
	--samples-file $samples \
	--regions-file $regions \
	--cromwell

