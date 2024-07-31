#!/bin/bash

chr="chr15"
ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz"
#ref="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/phased_ACAN_vntr_snp_650_samples.sorted.vcf.gz"

# VCF file with all samples
gt="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.${chr}.vcf.bgz"
samples="$WORKSPACE_BUCKET/saraj/vntr_samples/sample_ids.txt"
regions="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_50m.bed"
map="$WORKSPACE_BUCKET/saraj/genetic_map/plink.chr15_w_chr.GRCh38.map"

        #--name ${chr}_acan_50mbp_aou_50k_srwgs_samples \
time python imputation_aou.py \
        --vcf $gt \
        --ref-panel $ref \
        --name batch_set1 \
        --window 20 \
	--overlap 2 \
	--chrom $chr \
	--mem 80 \
	--map $map \
	--samples-file $samples \
	--regions-file $regions \
	--batch-size 1000
	#--dryrun \
	#--cromwell


