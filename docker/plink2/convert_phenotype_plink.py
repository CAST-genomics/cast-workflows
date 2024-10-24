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
import gcsfs


    # Get token and set up project
token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
    capture_output=True, check=True, encoding='utf-8')
token = str.strip(token_fetch_command.stdout)

#debug project"
project = os.getenv("GCS_REQUESTER_PAYS_PROJECT")

print(f"GOOGLE_PROJECT: {project}")


project_fetch_command = subprocess.run('echo $GOOGLE_PROJECT', shell=True, 
                                       capture_output=True, check=True, encoding='utf-8')
project_fetch = str.strip(project_fetch_command.stdout)
print(f"GOOGLE_PROJECT_fetch: {project_fetch}")

print("Current Environment Variables:")
print(os.environ)  # Print all environment variables



def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def DownloadPT(filename):
    """
	Download phenotype_covar.csvlocally

	Arguments
	---------
	filename : str
	   GCP path
	"""
    cmd = "gsutil cp {filename} .".format(filename=filename)
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    print(output.decode("utf-8"))
    return filename.split("/")[-1]  # Return local filename


def GetFloatFromPC(x):
    x = x.replace("[","").replace("]","")
    return float(x)

def LoadAncestry(ancestry_pred_path,token,project):
    fs = gcsfs.GCSFileSystem(token=token,project=project,requester_pays=True)
    with fs.open(ancestry_pred_path, 'r') as file:
        ancestry =  pd.read_csv(ancestry_pred_path, sep="\t")
    #if ancestry_pred_path.startswith("gs://"):
    #    if not os.path.isfile("ancestry_preds.tsv"):
    #        os.system("gsutil -u %s cp %s ."%(project,ancestry_pred_path))
    #    ancestry_pred_path = "ancestry_preds.tsv"
    #ancestry_pred_path = "ancestry_preds.tsv"
    ancestry.rename({"research_id": "IID"}, axis=1, inplace=True)
    ancestry['IID'] = ancestry['IID'].astype(str)
    num_pcs = len(ancestry["pca_features"].values[0].split(","))
    pcols = ["PC_%s"%i for i in range(num_pcs)]
    ancestry[pcols] = ancestry["pca_features"].str.split(",", expand=True)
    for p in pcols:
        ancestry[p] = ancestry[p].apply(lambda x: GetFloatFromPC(x), 1)
    ancestry.insert(0, 'FID', 0)

    return ancestry

def convert_csv_to_plink (ptfile):
    df = pd.read_csv(ptfile)
    df.insert(0, 'FID', 0)
    df.rename(columns={"person_id": "IID"}, inplace=True)
    return df

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions",type=str, default="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
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
    covars = pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    ancestry = LoadAncestry(args.ancestry_pred_path,token,project)
    
    print(ancestry)
    plink = convert_csv_to_plink(DownloadPT(ptcovar_path))
    print(plink)

    plink['IID'] = plink['IID'].astype(str)

    data = pd.merge(plink[["FID","IID"]+covars], ancestry[["IID"]+pcols],on=["IID"],how="inner")
    plink_pheno = plink[["FID","IID","phenotype"]]
    plink_pheno.to_csv(f"{args.phenotype}_pheno_plink.txt", sep="\t", index=False)
    data.to_csv(f"{args.phenotype}_covar_combined.txt", sep="\t", index=False)
    
    sys.exit(0)
    print(f"Done converting {args.phenotype} to plink format")


if __name__ == "__main__":
    main()


