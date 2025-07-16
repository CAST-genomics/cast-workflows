#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 
chrom=21
./run_convert_to_pgen.py \
--name SNP_imputation_chr${chrom} \
--chrom ${chrom} 

"""

# need to change GetFileBatches
import argparse
import json
import os
import subprocess
import sys

sys.path.append("../utils")
import aou_utils

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the output file", required=True, type=str)
	parser.add_argument("--chrom", help="Which chromosome to process", type=str, required=True)
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
	#pfile = "tr_imputation/enstr-v3/results-250K/reannotated/"

	# Set up workflow JSON
	json_dict = {}
	json_dict = {}
	json_dict["find_peaks.pgens"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.pgen')]
	json_dict["find_peaks.psams"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.psam')]
	json_dict["find_peaks.pvars"] = [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.pvar')]
	json_dict["find_peaks.out_prefix"] = args.name


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
	aou_utils.RunWorkflow("wdl/convert_vcf_to_pgen.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
