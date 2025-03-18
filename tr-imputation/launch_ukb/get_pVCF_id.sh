#!/usr/bin/env bash

## This section is used to extract pvcf_file_ids
path_prefix="/Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]"



for chrom in $(seq 1 22);
do
  dx ls -l "${path_prefix}/chr${chrom}/*gz" | sort -t "_" -k3,3V | sed -n 's/.*(\(.*\)).*/\1/p' > ./chr${chrom}_pvcf_file_ids.txt
done



# This section is trim vcf using the https://github.com/drarwood/vcf_trimmer/tree/master?tab=readme-ov-file

#chrom=21
#out_path="/yal084/tr_imputation/trimmed_vcfs"
#applet_name="applet-GzF40ZQJX3JzyJgXyj0j5bjJ"

#dx rm -r ${out_path}/chr${chrom}/
#dx mkdir -p ${out_path}/chr${chrom}/batched_vcf_list/
#echo ${path_prefix}/chr${chrom}/*gz
#dx ls -l "${path_prefix}/chr${chrom}/*gz" | sort -t "_" -k3,3V | awk -v prefix="${path_prefix}/chr${chrom}/" 'newpath=prefix$0 {print newpath}' > ./chr${chrom}_pvcf_list_full_path.txt
# ./chr${chrom}_pvcf_list.txt >./chr${chrom}_pvcf_list_full_path.txt

# mkdir ./chr${chrom}
# split -l 200 -d -a 4 ./chr${chrom}_pvcf_list_full_path.txt ./chr${chrom}/chr${chrom}_pvcf_batch_
#
# rm_gt_fields="FORMAT/FT,FORMAT/AD,FORMAT/MD,FORMAT/DP,FORMAT/RA,FORMAT/PP,FORMAT/GQ"
# qc_thresholds="INFO/AAScore>=0.5"
#
# cd ./chr${chrom}
# for batch_file in $(ls ./ | head -2); do
#   dx upload ${batch_file} --path ${out_path}/chr${chrom}/batched_vcf_list/
#
#   echo "submit jobs for $(basename ${batch_file})" 
#   dx run ${applet_name} \
#     -ivcf_file_list=${out_path}/chr${chrom}/batched_vcf_list/$(basename ${batch_file}) \
#     -ifile_label=trimmed \
#     -ioutput_dir=${out_path}/chr${chrom}/ \
#     -iqc_thresholds=${qc_thresholds} \
#     -ifields_to_remove=${rm_gt_fields} \
#     -iconcurrent_processes=10 \
#     -ithreads=2 \
#     -y
#   done
