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

sys.path.append("../utils")
import aou_utils

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="name of the run", required=True, type=str)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
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
	pfile = "tr_imputation/enstr-v3/results-250K/"

	
	# Set up workflow JSON
	json_dict = {}
	json_dict["tr_gwas.pgens"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.pgen')]
	json_dict["tr_gwas.psams"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.psam')]
	json_dict["tr_gwas.pvars"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.pvar')]
	json_dict["tr_gwas.phenotypes"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix="phenotypes/") if blob.name.endswith('.csv')]
	json_dict["tr_gwas.cohorts"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix="samples/") if blob.name.endswith('.csv')]
	
	

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
