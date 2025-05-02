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
    for bfile in ["bgen", "sample"]:
        curr_file = f"ukb22828_c{chrom}_b0_v3.{bfile}"
        results = list(dxpy.find_data_objects(name=f"{curr_file}", folder=f"{file_path}", project=f"{project_id}", return_handler=True))
        file_id.append(results[0].get_id())
    return file_id


def saveJson(file_name, input_dict):
    with open(file_name, "w") as f:
        json.dump(input_dict, f, indent=4)


def RunWorkflow(json_dict, workflow_id, chrom, job_priority="high"):
    workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
    analysis = workflow.run(json_dict, folder=f"/yal084/QCed_pgen_from_field_22828/chr{chrom}", priority=job_priority)
    sys.stderr.write(f"Submitted {analysis.get_id()} for chr{chrom}") 
    return analysis

# plink_options = "--mac 10 --maf 0.0001 --hwe 1e-15 --mind 0.1 --geno 0.1"

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--chrom", help="chromosome need to process, e.g. 21", required=True, type=int)
    parser.add_argument("--bfile-path", help="the folder of bfile on DNANexus", required=False, default="/Bulk/Imputation/UKB imputation from genotype", type=str)
    parser.add_argument("--split-mem", help="memory of instance used for plist", required=False, default=32, type=int)
    parser.add_argument("--bgen-qc", help="qc for bgen file, default from UKB recommended", required=False, default="--mac 20", type=str)
    parser.add_argument("--workflow-id", help="id of workflow", required=False, default="workflow-J0BKGx8JX3JjPYjQJ5F1Yvfv", type=str)

    args = parser.parse_args()
    bgen_id, sample_id = getbfileID(args.bfile_path, args.chrom)    

    json_dict = {}
    json_dict["stage-common.bgen_file"] = dxpy.dxlink(bgen_id)
    json_dict["stage-common.sample_file"] = dxpy.dxlink(sample_id)
    json_dict["stage-common.split_mem"] = args.split_mem
    json_dict["stage-common.chrom"] = f"chr{args.chrom}"
    json_dict["stage-common.bgen_qc"] = args.bgen_qc 

    saveJson(f"./logs/chr{args.chrom}_bgenToQCpgen.json", json_dict)
    RunWorkflow(json_dict, args.workflow_id, args.chrom)


if __name__ == "__main__":
    main()

