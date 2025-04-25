#!/usr/bin/env bash

## This section is used to extract pvcf_file_ids
path_prefix="/Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]"
genotype_batch_prefix="/yal084/vcf_from_genotype_bfile/"

mkdir -p ./genotype_batched_file_id/
for chrom in $(seq 1 22);
do
  echo ${chrom}
  dx ls -l "${genotype_batch_prefix}/chr${chrom}/liftOver_VCFs/splitted_by_samples/*gz" | sort -t "_" -k3,3V | sed -n 's/.*(\(.*\)).*/\1/p' > ./genotype_batched_file_id/chr${chrom}_batchVCF_file_ids.txt
  dx ls -l "${genotype_batch_prefix}/chr${chrom}/liftOver_VCFs/splitted_by_samples/*.tbi" | sort -t "_" -k3,3V | sed -n 's/.*(\(.*\)).*/\1/p' > ./genotype_batched_file_id/chr${chrom}_batchVCF_idx_ids.txt

done
