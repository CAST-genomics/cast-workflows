#!/usr/bin/env python3

"""
Preprocess AoU case control phenotypes.
Output cleaned phenotype and covariates files
for the specified phenotype

Requires these environment variables to be set:
WORKSPACE_BUCKET
"""

import argparse
import os
import pandas as pd
import subprocess
import sys

bucket = os.environ["WORKSPACE_BUCKET"]
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")

def MSG(msg_str):
    """
    Write a helpful progress message to stderr

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    sys.stderr.write("[aou_phenotype_preprocessing]: %s\n"%msg_str.strip())

def ERROR(msg_str):
    """
    Write an error message to stderr then quit
    with exit code 1

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    MSG("ERROR: " + msg_str)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phecodeX", help="phecodeX for phenotype", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs,sex to keep", type=str, default=SAMPLEFILE)

    args = parser.parse_args()
    MSG("Processing %s"%args.phecodeX)

    # Check if valid phecodeX
    subprocess.run("gsutil cp "+bucket+"/phenotypes/allbyall_phecodeX_headers.csv ./", shell=True, check=True)
    headers = pd.read_csv("allbyall_phecodeX_headers.csv")
    try:
        phecode_idx = headers[headers['phecodeX'] == args.phecodeX].index[0]
    except IndexError:
        ERROR(f"PhecodeX {args.phecodeX} is not supported by this script.")
    phecode_idx += 1 #awk index starts at 1

    # Get phenotype and covariates
    if not os.path.isfile("mcc2_phecodex_table.csv"):
        # This google bucket is for the AOU all by all PhecodeX workspace
        # https://workbench.researchallofus.org/workspaces/aou-rw-2e6e1444/allbyallphecodexphenotypescuration/data
        phecode_table_loc = "gs://fc-secure-2607142e-b58a-40b2-9174-5cfed82a428b/notebooks/phenotype_data/mcc2_phecodex_table.csv"
        subprocess.run(f"gsutil cp {phecode_table_loc} ./", shell=True, check=True)
    awk_command = f"awk -F',' '{{print $1 \",\" ${phecode_idx}}}' mcc2_phecodex_table.csv > {args.phecodeX}.csv"
    subprocess.run(awk_command, shell=True, check=True)
    ptdata = pd.read_csv(f"{args.phecodeX}.csv")
    ptdata = ptdata[ptdata[args.phecodeX].notnull()]
    ptdata[args.phecodeX] = ptdata[args.phecodeX].astype('int')

    subprocess.run("gsutil cp "+bucket+"/phenotypes/allbyall_demographics_table.csv ./", shell=True, check=True)
    demog = pd.read_csv("allbyall_demographics_table.csv")
    demog = demog[["person_id","age_at_cdr"]]
    demog.rename(columns={"age_at_cdr":"age"}, inplace=True)
    

    data = pd.merge(ptdata, demog, on="person_id", how="inner")
    MSG("After merge, have %s data points"%data.shape[0])

    # Restrict to samples we want to keep
    sampfile = args.samples
    if sampfile.startswith("gs://"):
        sampfile = sampfile.split("/")[-1]
        if not os.path.isfile(sampfile):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(args.samples))
    samples = pd.read_csv(sampfile)
    data = data.merge(samples, on='person_id', how='inner')
    MSG("After filter samples, have %s data points"%data.shape[0])

    # Output final phenotype value
    data.rename(columns={args.phecodeX: "phenotype"}, inplace=True)
    data[["person_id", "phenotype", "age", "sex_at_birth_Male"]].to_csv(args.phecodeX+"_phenocovar.csv", index=False)

if __name__ == "__main__":
    main()
