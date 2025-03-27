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

def getfileID(file_path, chrom, project_id="project-GbjB4x0JX3JyK5PZK67q0fZG"):
    
    map_file =  f"plink.chr{chrom}.GRCh38_new.map"
    result =  list(dxpy.find_data_objects(name=f"{map_file}", folder=f"{file_path}/genetic_maps/", project=f"{project_id}", return_handler=True))
    map_id = result[0].get_id()
    ref_file = f"ensembletr_refpanel_v4_chr{chrom}.bref3"
    result =  list(dxpy.find_data_objects(name=f"{ref_file}", folder=f"{file_path}/ensembletr_refpanel_v4/", project=f"{project_id}", return_handler=True))
    ref_id = result[0].get_id()
    vcf_file = f"ensembletr_refpanel_v4_chr{chrom}.vcf.gz"
    result =  list(dxpy.find_data_objects(name=f"{vcf_file}", folder=f"{file_path}/ensembletr_refpanel_v4/", project=f"{project_id}", return_handler=True))
    vcf_id = result[0].get_id()
    index_file = f"ensembletr_refpanel_v4_chr{chrom}.vcf.gz.tbi"
    result =  list(dxpy.find_data_objects(name=f"{index_file}", folder=f"{file_path}/ensembletr_refpanel_v4/", project=f"{project_id}", return_handler=True))
    index_id = result[0].get_id()
    return map_id, ref_id, vcf_id, index_id 


def saveJson(file_name, input_dict):
    with open(file_name, "w") as file:
        json.dump(input_dict, file, indent=4)

def RunWorkflow(json_dict, workflow_id, name, depends=[], run_priority="normal"):
    
    workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
    analysis = workflow.run(json_dict,
                            depends_on=[item.get_id() for item in depends],
                            folder=f"/yal084/imputed_files/{name}",
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
    parser.add_argument("--chrom", help="chromosome of pVCF to trim, e.g. 1", required=True, type=int)
    parser.add_argument("--imputation-mem", help="memory for beagle imputation", required=False, default=32, type=int) 
    parser.add_argument("--imputation-id", help="workflow-id for beagle imputation", required=False, default="workflow-GzX4fz8JX3Jg5fV9vp6KPJ0z", type=str)
    parser.add_argument("--file-path", help="the folder of bfile on DNANexus", required=False, default="/additional_files", type=str)
    parser.add_argument("--batch-size", help="number of VCF files run in a job", default=10, type=int, required=False)
    parser.add_argument("--batch-num", help="max number of jobs submitted, for test", required=False, default=-1, type=int)
    
    parser.add_argument("--extract-id", help="workflow for extract SNP/TRs", required=False, default="workflow-GzX83B0JX3JqxvBb4vJJ0j5b", type=str)
    parser.add_argument("--extract-mem", help="memory for extraction", required=False, default=16, type=int)
    parser.add_argument("--extract-merge-mem", help="memory for merge", required=False, default=32, type=int)

    parser.add_argument("--merge-anno-id", help="workflow for merge and annotation", required=False, default="workflow-GzX96z8JX3JgFV97PXp8Jq0q", type=str)
    parser.add_argument("--merge-mem", help="memory for extraction", required=False, default=32, type=int)
    parser.add_argument("--anno-mem", help="memory for merge", required=False, default=32, type=int)

    args = parser.parse_args()
    log_dir="./logs/batch_imputation/"
    os.makedirs(log_dir, exist_ok=True)
    map_id, ref_bref3_id, ref_vcf_id, ref_index_id = getfileID(args.file_path, args.chrom)

    impute_dict = {}
    impute_dict["stage-common.ref_panel"] = dxpy.dxlink(ref_bref3_id)
    impute_dict["stage-common.genetic_map"] = dxpy.dxlink(map_id) 
    impute_dict["stage-common.mem"] = args.imputation_mem 
    
    sys.stderr.write("Setting up batches...\n")
    snp_batches = GetFileBatches(f"./genotype_batched_file_id/chr{args.chrom}_batchVCF_file_ids.txt", args.batch_size, args.batch_num)
    index_batches = GetFileBatches(f"./genotype_batched_file_id/chr{args.chrom}_batchVCF_idx_ids.txt", args.batch_size, args.batch_num)   
    
    depends = []
    curr_idx = 0
    while curr_idx < len(snp_batches):
        batch_dict = impute_dict.copy()
        batch_dict["stage-common.snp_vcf_list"] = snp_batches[curr_idx]
        batch_dict["stage-common.snp_idx_list"] = index_batches[curr_idx]
        impute_job = RunWorkflow(batch_dict, args.imputation_id, f"impute_from_genotyped/raw_imputed_files/chr{args.chrom}/", depends=[]) # , run_priority="normal")
        saveJson(f"{log_dir}/imputation_chr{args.chrom}_CHUNK{curr_idx}.json", batch_dict)
        depends.append(impute_job)
        curr_idx += 1

    # extract TRs
    sys.stderr.write(f"Setting up extraction ...\n")
    extract_job_list = [] 
    for impute_job, curr_idx in zip(depends, range(0, curr_idx)):
        vcf = GetJBOR(impute_job, "stage-outputs.imputed_snp_strs")
        vcf_idx = GetJBOR(impute_job, "stage-outputs.imputed_snp_strs_idx")
        # print(vcf_idx)

        extract_dict = {}
        extract_dict["stage-common.vcf_list"] = vcf 
        extract_dict["stage-common.idx_list"] = vcf_idx 
        extract_dict["stage-common.extract_mem"] = args.extract_mem
        extract_dict["stage-common.merge_mem"] = args.extract_merge_mem
        extract_dict["stage-common.merge_prefix"] = f"chr{args.chrom}_CHUNK{curr_idx}"

        extract_job = RunWorkflow(extract_dict, args.extract_id, f"impute_from_genotyped/TR_merged_by_chunks/chr{args.chrom}/", depends=depends) #, run_priority="high")
        extract_job_list.append(extract_job) 
        saveJson(f"{log_dir}/extract_chr{args.chrom}_CHUNK{curr_idx}.json", extract_dict)
   
    # setup merge and annotation
    chunk_tr_vcfs = []
    chunk_tr_vcfs_idx = []
    sys.stderr.write(f"Setting up merge and annotation...\n")
    for extract_job in extract_job_list:
        vcf = GetJBOR(extract_job, "stage-outputs.vcf_files")
        vcf_idx = GetJBOR(extract_job, "stage-outputs.index_files")
        chunk_tr_vcfs.append(vcf)
        chunk_tr_vcfs_idx.append(vcf_idx)

    anno_dict = {}
    anno_dict["stage-common.vcf_list"] = chunk_tr_vcfs
    anno_dict["stage-common.idx_list"] = chunk_tr_vcfs_idx
    anno_dict["stage-common.ref_vcf"] = dxpy.dxlink(ref_vcf_id)
    anno_dict["stage-common.ref_vcf_idx"] = dxpy.dxlink(ref_index_id)
    anno_dict["stage-common.merge_mem"] = args.merge_mem 
    anno_dict["stage-common.out_prefix"] = f"chr{args.chrom}" 
    anno_dict["stage-common.anno_mem"] = args.anno_mem
    RunWorkflow(anno_dict, args.merge_anno_id, f"impute_from_genotyped/annotated_TRs/chr{args.chrom}/", depends=extract_job_list)
    saveJson(f"{log_dir}/merge_and_anno_chr{args.chrom}.json", anno_dict)

if __name__ == "__main__":
    main()

