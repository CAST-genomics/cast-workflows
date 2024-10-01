#!/usr/bin/env python3



"""
Convert AoU phenotypes to plink friendly format.
Output cleaned phenotype and comined covariates files
for the specified phenotype

Usage:
./convert_phenotype_plink.py --phenotype <phenotype>

"""


import os
import pandas as pd
import subprocess
import sys
import csv
import argparse
import os


ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"


def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def GetFloatFromPC(x):
    x = x.replace("[","").replace("]","")
    return float(x)

def LoadAncestry(ancestry_pred_path):
    if ancestry_pred_path.startswith("gs://"):
        if not os.path.isfile("ancestry_preds.tsv"):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(ancestry_pred_path))
        ancestry_pred_path = "ancestry_preds.tsv"
    ancestry = pd.read_csv(ancestry_pred_path, sep="\t")
    ancestry.rename({"research_id": "IID"}, axis=1, inplace=True)
    ancestry['IID'] = ancestry['IID'].astype(str)
    num_pcs = len(ancestry["pca_features"].values[0].split(","))
    pcols = ["PC_%s"%i for i in range(num_pcs)]
    ancestry[pcols] = ancestry["pca_features"].str.split(",", expand=True)
    for p in pcols:
        ancestry[p] = ancestry[p].apply(lambda x: GetFloatFromPC(x), 1)
    ancestry.insert(0, 'FID', 0)

    return ancestry

def convert_csv_to_plink (ptfile,outfile):
    with open(ptfile, 'r') as csv_file, open(outfile, 'w', newline='') as plink_file:
        reader = csv.reader(csv_file)
        writer = csv.writer(plink_file, delimiter='\t')
        
        for i, row in enumerate(reader):
            if i == 0:
                # Print header with "FID", "IID", and the second column
                writer.writerow(["FID", "IID"], row[1:])
            else:
                # Print "0", first column, and second column
                writer.writerow(["0", row[0], row[1], row[2], row[3]])
    
    return outfile

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")

    args = parser.parse_args()

    # Set up paths
    if args.phenotype.endswith(".csv"):
        ptcovar_path = args.phenotype
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)


    # Get covarlist
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""]
    pt_covars = [item for item in args.ptcovars.split(",") if item != ""]
    covars = pcols + pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    data = pd.read_csv(ptcovar_path)
    ancestry = LoadAncestry(ANCESTRY_PRED_PATH)
    pheno_file = convert_csv_to_plink (data,f"{args.phenotype}_plink")
    plink_pheno = pd.read_csv(pheno_file,sep='\t')
    data = pd.merge(plink_pheno, ancestry[["IID"]+covars], on=["IID"],how="inner")

    data.to_csv(f"{args.phenotype}_covar_combined.txt", sep="\t", index=False)



if __name__ == "__main__":
    main()


