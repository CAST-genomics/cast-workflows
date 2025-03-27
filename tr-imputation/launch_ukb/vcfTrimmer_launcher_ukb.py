#!/usr/bin/env python3

import argparse
import dxpy
import sys
import time
import json
import os

def GetFileBatches(file_list, batch_size, batch_num=-1):
    """


    """
    vcf_batches = []
    curr_batch_vcfs = []

    with open(file_list, "r") as f:
        for line in f:
            if len(curr_batch_vcfs) == batch_size:
                vcf_batches.append(curr_batch_vcfs)
                if len(vcf_batches) == batch_num: break
                curr_batch_vcfs = []
            vcf_id = line.strip()
            curr_batch_vcfs.append(dxpy.dxlink(vcf_id))

    if len(curr_batch_vcfs) > 0 and (len(vcf_batches) < batch_num or batch_num == -1):
        vcf_batches.append(curr_batch_vcfs)
    return vcf_batches 


def saveJson(file_name, input_dict):
    with open(file_name, "w") as file:
        json.dump(input_dict, file, indent=4)

def RunWorkflow(json_dict, workflow_id, name, depends=[], run_priority="low"):
    
    workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
    analysis = workflow.run(json_dict,
                            depends_on=[item.get_id() for item in depends],
                            folder=f"/yal084/trimmed_pvcfs/{name}",
                            priority=run_priority)
    sys.stderr.write(f"Started analysis {analysis.get_id()} job: {name}\n")
    time.sleep(30)
    return analysis


def GetJBOR(analysis, filename):
    return {
        "$dnanexus_link": {
            "analysis": analysis.get_id(),
            "field": filename
        }
    }

def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--chrom", help="chromosome of pVCF to trim, e.g. chr1", required=True, type=str)
    parser.add_argument("--batch-size", help="number of pVCF files run in a job", default=10, type=int, required=False)
    parser.add_argument("--batch-num", help="max number of jobs submitted, for test", required=False, default=-1, type=int)
    parser.add_argument("--qc-thresholds", help="QC for Graphertype calls", required=False, default="INFO/AAScore>=0.5", type=str)
    parser.add_argument("--rm-tags", help="tags to remove/keep, follow bcftools format", required=False, default="^INFO/AAScore,INFO/AC,INFO/AF,^FORMAT/GT, FORMAT/PL", type=str)
    parser.add_argument("--threads-num", help="no. of threads for bcftools", required=False, default=2, type=int)
    parser.add_argument("--trimvcf-mem", help="memory for trim vcf", required=False, default=16, type=int) 
    parser.add_argument("--trimvcf-id", help="workflow-id to trim vcfs", required=False, default="workflow-GzVZz8jJX3Jyk8jjzkpBgyqK", type=str)

#    parser.add_argument("--concat-id", help="workflow-id to concat vcfs", required=False, default="workflow-GzVY5B0JX3JpX6bB8JPJqP5g", type=str)
    parser.add_argument("--concat-id", help="workflow-id to concat vcfs", required=False, default="workflow-GzXFpXjJX3Jb334BVYgZGygg", type=str)

    parser.add_argument("--concat-mem", help="memory for concat vcfs", required=False, default=64, type=int)

    parser.add_argument("--splitvcf-id", help="workflow-id to split vcfs", required=False, default="workflow-GzQybp0JX3JQjX2j6VzfbzKx", type=str)
    parser.add_argument("--splitvcf-mem", help="memory for split vcfs", required=False, default=64, type=int)
    parser.add_argument("--sample-size", help="number of samples in each batch", required=False, default=1000, type=int)

    args = parser.parse_args()
    log_dir="./logs/trimvcfs/"
    os.makedirs(log_dir, exist_ok=True)
    trim_dict = {}
    # trim_dict["stage-common.outprefix"] = args.chrom
    trim_dict["stage-common.qc_thresholds"] = args.qc_thresholds
    trim_dict["stage-common.rm_tags"] = args.rm_tags
    trim_dict["stage-common.threads_num"] = args.threads_num
    trim_dict["stage-common.bcftools_mem"] = args.trimvcf_mem

    sys.stderr.write("Setting up batches...\n")
    vcf_batches = GetFileBatches(f"./pvcf_id/{args.chrom}_pvcf_file_ids.txt", args.batch_size, args.batch_num)
    
    depends = []
    curr_idx = 0
    while curr_idx < len(vcf_batches):
        batch_name = f"{args.chrom}-CHUNK{curr_idx}"
        batch_dict = trim_dict.copy()
        batch_dict["stage-common.outprefix"] = batch_name
        batch_dict["stage-common.pvcf_batches"] = vcf_batches[curr_idx]
        trim_job = RunWorkflow(batch_dict, args.trimvcf_id, f"{args.chrom}/by_CHUNKs/", depends=[], run_priority="high")
        saveJson(f"{log_dir}/trimVCF_{args.chrom}_CHUNK{curr_idx}.json", batch_dict)
        depends.append(trim_job)
        curr_idx += 1


    concat_vcfs = []
    concat_vcfs_idx = []
    for trim_job in depends:
        vcf = GetJBOR(trim_job, "stage-outputs.final_vcf")
        vcf_idx = GetJBOR(trim_job, "stage-outputs.final_vcf_indx")
        concat_vcfs.append(vcf)
        concat_vcfs_idx.append(vcf_idx)

    sys.stderr.write("Setting up concat from all-batches...\n")
    concat_dict = {}
    concat_dict["stage-common.merged_prefix"] = args.chrom
    concat_dict["stage-common.vcf_files"] = concat_vcfs
    concat_dict["stage-common.vcf_indexs"] = concat_vcfs_idx
    concat_dict["stage-common.mem"] = args.concat_mem
    
    concat_job = RunWorkflow(concat_dict, args.concat_id, args.chrom, depends=depends, run_priority="high")
    saveJson(f"{log_dir}/concat_{args.chrom}.json", concat_dict)


    split_dict = {}
    combined_vcf = GetJBOR(concat_job, "stage-outputs.final_vcf")
    split_dict["stage-common.vcf_file"] = combined_vcf
    split_dict["stage-common.chromosome"] = args.chrom
    split_dict["stage-common.sample_size"] = args.sample_size
    split_dict["stage-common.split_mem"] = args.splitvcf_mem

    RunWorkflow(split_dict, args.splitvcf_id, f"{args.chrom}/splitted_by_samples/", depends=[concat_job], run_priority="high")
    saveJson(f"{log_dir}/splitVCF_{args.chrom}.json", split_dict)

if __name__ == "__main__":
    main()

