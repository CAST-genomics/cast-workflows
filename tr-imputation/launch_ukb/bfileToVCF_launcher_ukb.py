#!/usr/bin/env python3

'''
Scripts that convert bfile to vcf for imputation
'''

import argparse
import dxpy
import sys
import json

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


def RunWorkflow(json_dict, workflow_id, chrom):
    workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
    analysis = workflow.run(json_dict, folder=f"/yal084/vcf_from_genotype_bfile/chr{chrom}")
    sys.stderr.write(f"Submitted {analysis.get_id()} for chr{chrom}") 
    return analysis

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--chrom", help="chromosome need to process, e.g. 21", required=True, type=int)
    parser.add_argument("--bfile-path", help="the folder of bfile on DNANexus", required=False, default="/Bulk/Genotype Results/Genotype calls", type=str)
    parser.add_argument("--ref-genome", help="reference genome for check REF/ALT in bfile", required=False, default="file-GzKvKp8JX3JvVjg9Qp4pzqgv", type=str)
    parser.add_argument("--batch-size", help="number of samples per batch", required=False, default=1000, type=int)
    parser.add_argument("--split-mem", help="memory of instance used for plist", required=False, default=64, type=int)
    parser.add_argument("--workflow-id", help="id of workflow", required=False, default="workflow-GzKy00jJX3JyY8PqQQXff99Y", type=str)
 
    args = parser.parse_args()
    bim_id, bed_id, fam_id = getbfileID(args.bfile_path, args.chrom)    

    json_dict = {}
    json_dict["stage-common.bed"] = dxpy.dxlink(bed_id)
    json_dict["stage-common.bim"] = dxpy.dxlink(bim_id)
    json_dict["stage-common.fam"] = dxpy.dxlink(fam_id)
    json_dict["stage-common.ref_genome"] = dxpy.dxlink(args.ref_genome)
    json_dict["stage-common.batch_size"] = args.batch_size
    json_dict["stage-common.split_mem"] = args.split_mem
    json_dict["stage-common.chrom"] = f"chr{args.chrom}"
    saveJson(f"./logs/chr{args.chrom}_bfileToVCF.json", json_dict)
    RunWorkflow(json_dict, args.workflow_id, args.chrom)


if __name__ == "__main__":
    main()

