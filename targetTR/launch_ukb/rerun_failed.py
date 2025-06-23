#!/usr/bin/env python3

import os
import dxpy
import time
import json
import sys
import argparse
import dxpy
from datetime import date
from glob import glob


def RunWorkflow(json_dict, workflow_id, name, out_dir, job_priority, depends=[]):
    """
    Run DNA Nexus workflow

    Arguments
    ---------
    json_dict : dict
        Input arguments for the workflow
    workflow_id : str
        ID of the DNA Nexus worfklow
    name : str
        Used to determine where to store output

    Returns
    -------
    analysis : DXAnalysis object
    """
    workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
    analysis = workflow.run(json_dict, \
                            depends_on=[item.get_id() for item in depends], \
                            folder=f"{out_dir}", \
                            name=name.split("/")[-1], 
                            priority=job_priority)
    sys.stderr.write("Started analysis %s (%s)\n"%(analysis.get_id(), name))
    # sleep to delay launches between analyses as each analysis mounts
    # a bunch of files which can put pressure on the API call limits
    # imposed by DNANexus that they care about but don't cleanly enforce
    time.sleep(30)
    return analysis


def loadJson(json_file):
    with open(json_file, "r") as f:
        job_dict = json.load(f)
        return job_dict


def saveJason(file, input_dict):
    with open(file, "w") as f:
        json.dump(input_dict, f, indent=4)


def main():
    parser = argparse.ArgumentParser(__doc__)
    #parser.add_argument("--rerun-file", help="path to file containing list of re-run jobs", required=True, type=str)
    parser.add_argument("--job-priority", help="use spot instance or not", required=False, default="normal", type=str)
    parser.add_argument("--save-path", help="fold to save on DNANexus", default="/TargetedSTR/results", required=False, type=str)
    parser.add_argument("--json-path", help="path to access json file for job launching", required=True, type=str)
    parser.add_argument("--targetTR-workflow", help="workflow id of targetTR", required=False, default="workflow-GpkxfqjJv7B4pGXPp5K9V6q2", type=str)
    args = parser.parse_args()
    today = date.today()
    formatted_date = today.strftime("%y%m%d")
    log_dir = f"./submission_logs/{formatted_date}_rerun/"
    print(f"rerun logs stored at {log_dir}")
    os.makedirs(log_dir, exist_ok=True)
    workflow_id = args.targetTR_workflow
    fail_num = 0
    try:
        job_id_file = glob(f"{args.json_path}/*analysis_id.txt")[0]
        with open(job_id_file, "r") as f:
            for line in f:
                analysis_id = line.strip()
                #print(analysis_id)
                analysis = dxpy.DXAnalysis(analysis_id)
                desc = analysis.describe()
                if desc["state"] != "done":
                    fail_num += 1
                    print(f"{fail_num:,}: {analysis_id} is failed")
                    job_name = desc["name"]
                    json_file = f"{args.json_path}/{job_name}.json"
                    job_dict = loadJson(json_file)
                    new_analysis = RunWorkflow(job_dict, workflow_id, f"rerun_{job_name}", out_dir=f"{args.save_path}/{args.json_path.split("/")[-1]}/batched/{job_name}", job_priority=f"{args.job_priority}")
                    saveJason(f"{log_dir}/{job_name}.json", job_dict)
                    new_analysis.describe()
    except:
        print("There is analysis ID file!!")

if __name__ == "__main__":
    main()
