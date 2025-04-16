#!/usr/bin/env python3

'''
Scripts that convert bfile to vcf for imputation
    1. First convert the bfile to VCF file
    2. Use pyliftOver to lift the coordinate from hg19 to hg38
    3. split into small files, every 1_000 samples for imputation
'''

import argparse
import dxpy
import sys
import json
import os

def getbfileID(file_path, chrom, project_id="project-GbjB4x0JX3JyK5PZK67q0fZG"):
    file_id = []
    for bfile in ["bim", "bed", "fam"]:
        curr_file = f"ukb22418_c{chrom}_b0_v2.{bfile}"
        results = list(dxpy.find_data_objects(name=f"{curr_file}", folder=f"{file_path}", project=f"{project_id}", return_handler=True))
        file_id.append(results[0].get_id())
    return file_id


def saveJson(file_name, input_dict):
    with open(file_name, "w") as f:
        json.dump(input_dict, f, indent=4)


def GetJBOR(analysis, filename):
    return {
        "$dnanexus_link": {
            "analysis": analysis.get_id(),
            "field": filename
        }
    }


def RunWorkflow(json_dict, workflow_id, chrom, depends=[], job_priority="high"):
    workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
    analysis = workflow.run(json_dict, 
                            depends_on=[item.get_id() for item in depends],
                            folder=f"/yal084/vcf_from_genotype_bfile/chr{chrom}",
                            priority=job_priority)
    sys.stderr.write(f"Submitted {analysis.get_id()} for chr{chrom}...\n") 
    return analysis

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--bfile-path", help="the folder of bfile on DNANexus", required=False, default="/Bulk/Genotype Results/Genotype calls", type=str)
    parser.add_argument("--ref-genome", help="reference genome for check REF/ALT in bfile", required=False, default="file-GzKvKp8JX3JvVjg9Qp4pzqgv", type=str)
    parser.add_argument("--convert-mem", help="memory for converting the bfile to vcf", required=False, default=16, type=int)
    parser.add_argument("--chrom", help="chromosome need to process, e.g. 21", required=True, type=int)
    parser.add_argument("--bfileTovcf-id", help="id of workflow", required=False, default="workflow-GzV23bQJX3JYQVgKPBBqB943", type=str)
   
    # input for liftOver
    parser.add_argument("--liftover-mem", help="memory for liftOver", required=False, default=32, type=int)
    parser.add_argument("--liftover-id", help="workflow id for liftOvr", required=False, default="workflow-GzVQ9Q0JX3Jq44vFPpkGV4j8", type=str)
    
    # input for split by samples
    parser.add_argument("--split-mem", help="memory of instance used for plist", required=False, default=32, type=int)
    parser.add_argument("--splitVCF-id", help="workflow id for liftOvr", required=False, default="workflow-GzpbYx0JX3JXF743Jj3pbfgz", type=str)
    parser.add_argument("--sample-size", help="number of samples per batch", required=False, default=1000, type=int)
    parser.add_argument("--batch-num", help="number of scatter jobs", required=False, default=30, type=int)

    args = parser.parse_args()
    log_dir = "./logs/bfileToVCF"
    os.makedirs(log_dir, exist_ok=True)
    bim_id, bed_id, fam_id = getbfileID(args.bfile_path, args.chrom)    

    json_dict = {}
    json_dict["stage-common.bed"] = dxpy.dxlink(bed_id)
    json_dict["stage-common.bim"] = dxpy.dxlink(bim_id)
    json_dict["stage-common.fam"] = dxpy.dxlink(fam_id)
    json_dict["stage-common.ref_genome"] = dxpy.dxlink(args.ref_genome)
    json_dict["stage-common.convert_mem"] = args.convert_mem
    json_dict["stage-common.chrom"] = f"chr{args.chrom}"
    # json_dict["priority"] = "high" 

    bfileToVCF = RunWorkflow(json_dict, args.bfileTovcf_id, args.chrom, depends=[])
    saveJson(f"{log_dir}/chr{args.chrom}_bfileToVCF.json", json_dict)
    
    hg19_vcf = GetJBOR(bfileToVCF, "stage-outputs.converted_vcf")
    liftover_dict = {}
    liftover_dict["stage-common.hg19_vcf"] = hg19_vcf
    liftover_dict["stage-common.liftover_mem"] = args.liftover_mem
    # liftover_dict["priority"] = "high" 
    liftover = RunWorkflow(liftover_dict, args.liftover_id, f"{args.chrom}/liftOver_VCFs/", depends=[bfileToVCF]) 
    saveJson(f"{log_dir}/chr{args.chrom}_hg19Tohg38.json", liftover_dict)

    hg38_vcf = GetJBOR(liftover, "stage-outputs.hg38_vcf")
    split_dict = {}
    split_dict["stage-common.vcf_file"] = hg38_vcf
    split_dict["stage-common.chromosome"] = f"chr{args.chrom}"
    split_dict["stage-common.sample_size"] = args.sample_size 
    split_dict["stage-common.split_mem"] = args.split_mem
    split_dict["stage-common.batch_num"] = args.batch_num 

    # split_dict["priority"] = "high"
    RunWorkflow(split_dict, args.splitVCF_id, f"{args.chrom}/liftOver_VCFs/splitted_by_samples/", depends=[liftover], job_priority="normal")
    saveJson(f"{log_dir}/chr{args.chrom}_splitVCF.json", split_dict)

if __name__ == "__main__":
    main()
