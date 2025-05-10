#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 

./tr_gwas_aou.py 
--name 

"""

import argparse
import json
import os
import subprocess
import sys
from google.cloud import storage
import pandas as pd

sys.path.append("../utils")
import aou_utils

def DownloadAncestry(filename):
	"""
	Download a GCP path locally

	Arguments
	---------
	filename : str
	   GCP path
	"""
	cmd = "gsutil -u $GOOGLE_PROJECT cp {filename} .".format(filename=filename)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))

def GetPhenotypePath(phenotype, pheno_prefix, binary=False):
	"""
	Download a GCP path locally
	Arguments
	---------
	phenotype : str
	   GCP path
	binary: bool
	   If True, download case control phenotype
	"""
	phenotype_array = []
	phenotypes = [item.strip() for item in phenotype.split(',')]
	for item in phenotypes: 
		if binary:
			path = os.getenv("WORKSPACE_BUCKET")+"/phenotypes/case/"+pheno_prefix+item.strip()+"_phenocovar.csv"
			phenotype_array.append(path)
		else:
			path = os.getenv("WORKSPACE_BUCKET")+"/phenotypes/"+pheno_prefix+item.strip()+"_phenocovar.csv"
			phenotype_array.append(path)	
	return phenotype_array


def GetCohortPath(cohort):
	"""	
	Download a GCP path locally
	Arguments
	---------
	cohort : str		
	   GCP path
	"""
	cohort_array = []
	# Replace 'ALL' with 'passing_samples_v7.1' in cohort input and others
	cohort = cohort.replace("ALL", "passing_samples_v7.1")
	cohort = cohort.replace("AFR","AFR_BLACK")
	cohort = cohort.replace("EUR","EUR_WHITE")
	cohort = cohort.replace("AMR","AMR_HISPANIC")
	cohort = cohort.replace("NOT_AFR","NOT_AFR_BLACK")
	

	cohorts = [item.strip() for item in cohort.split(',')]
	for item in cohorts: 
		path = os.getenv("WORKSPACE_BUCKET")+"/samples/"+item.strip()+"_plink.txt"
		cohort_array.append(path)
	return cohort_array

def phenotypes_from_file(filename):
	pheno_df = pd.read_csv(filename)
	phenotypes = list(pheno_df["phecode_string"])
	return phenotypes


def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="name of the run", required=True, type=str)
	parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions",type=str, default="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--pheno-file", help="The phenotype file to extract phenotypes from.", type=str)
	parser.add_argument("--variant-as-covar-file", help="The file with variants that " + \
							"we want to include as covariate in association testing.", type=str)
	parser.add_argument("--phenotype", help="name of the phenotype, seperated by comma", required=False, type=str)
	parser.add_argument("--pheno-prefix", help="prefix of the phenocovar files if present", required=False, default="")
	parser.add_argument("--cohort", help="name of the cohort, seperated by comma, options: AFR, EUR, AMR, NOT_AFR, ALL", required=False, type=str)
	parser.add_argument("--logistic", help="Run logistic regression", action="store_true")
	args = parser.parse_args()


	# Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)
	#set up bucket
	bucket_name = os.getenv("WORKSPACE_BUCKET").replace("gs://", "")
	project = os.getenv("GOOGLE_PROJECT")
	client = storage.Client()
	bucket = client.bucket(bucket_name)
	# Define the gs prefix
	gs_prefix = f"gs://{bucket_name}/"
	#pfile = "tr_imputation/enstr-v3/results-250K/"
	pfile = "saraj/imputation_output/"

	variant_as_covar_file = "empty.csv"
	if args.variant_as_covar_file:
		variant_as_covar_file = args.variant_as_covar_file
	else:
		# Create an empty file if no file indicated.
		open(variant_as_covar_file, 'w').close()

	#return phenotype array if choose targeted phenotypes 
	if args.phenotype is not None:
			target_phenotype =  GetPhenotypePath(args.phenotype,
												binary=args.logistic,
												pheno_prefix=args.pheno_prefix)
	elif args.pheno_file:
		phenotypes = phenotypes_from_file(args.pheno_file)
		#print("len phenotypes extracted: ", len(set(phenotypes)))
		target_phenotype = []
		for phenotype in phenotypes:
			target_phenotype.extend(GetPhenotypePath(phenotype,
													 binary=args.logistic,
													 pheno_prefix=args.pheno_prefix))
	else:
		if args.logistic:
			target_phenotype = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix="phenotypes/case/") if blob.name.endswith('.csv')]
		else:
			target_phenotype = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix="phenotypes/") if blob.name.endswith('.csv')]
		if len(args.pheno_prefix) > 0:
			target_phenotype = [pheno for pheno in target_phenotype if args.pheno_prefix in pheno]

			

	print("len of target_phenotype: ", len(target_phenotype))
	print("pfile: ", pfile)
	#return cohort array if choose targeted cohorts 
	if args.cohort is not None:
		target_cohort =  GetCohortPath(args.cohort)
	else: 
		target_cohort = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix="samples/") if blob.name.endswith('.txt')]
	
	# Set up workflow JSON
	json_dict = {}
	json_dict["tr_gwas.pgens"] = [gs_prefix + blob.name for blob in client.list_blobs(bucket, prefix=pfile) if blob.name.endswith('.pgen')]
	json_dict["tr_gwas.psams"] = [gs_prefix + blob.name for blob in client.list_blobs(bucket, prefix=pfile) if blob.name.endswith('.psam')]
	json_dict["tr_gwas.pvars"] = [gs_prefix + blob.name for blob in client.list_blobs(bucket, prefix=pfile) if blob.name.endswith('.pvar')]
	json_dict["tr_gwas.cohorts"] = target_cohort
	json_dict["tr_gwas.GOOGLE_PROJECT"] = project
	json_dict["tr_gwas.GCS_OAUTH_TOKEN"] = token
	json_dict["tr_gwas.WORKSPACE_BUCKET"] = gs_prefix
	json_dict["tr_gwas.phenotypes"] = target_phenotype
	json_dict["tr_gwas.logistic"] = args.logistic
	json_dict["tr_gwas.variant_as_covar_file"] = variant_as_covar_file
	

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	aou_utils.RunWorkflow("tr_gwas.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)



if __name__ == "__main__":
	main()
