#!/usr/bin/env python3
"""
Script to launch 
 
./run_extract_SNP.py \
--name chr8 
--ID 6262ecbb-543a-4d2e-aedd-7f0c872ef3bc

"""

# need to change GetFileBatches
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
	parser.add_argument("--name", help="Name of the chromsome outfile", required=True, type=str)
	parser.add_argument("--ID", help="Name of the cromshell ID", type=str, required=True)
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
	vcfs = f"cromwell-execution/batch_imputation/{args.ID}/call-beagle/"

	# Set up workflow JSON
	json_dict = {}
	json_dict["extract_imputation_SNP.vcf"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=vcfs) if blob.name.endswith('.vcf.gz')]
	json_dict["extract_imputation_SNP.vcf_index"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=vcfs) if blob.name.endswith('.vcf.gz.tbi')]
	json_dict["extract_imputation_SNP.out_prefix"] = args.name


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
	aou_utils.RunWorkflow("wdl/extract_imputation_SNP.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
