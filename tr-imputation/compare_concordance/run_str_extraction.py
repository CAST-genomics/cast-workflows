#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 

./run_str_extraction.py \
--name ukb_imputation \
--str ukb_finemapped_hg38_str.txt


"""

# need to change GetFileBatches
import argparse
import json
import os
import subprocess
import sys
from google.cloud import storage

sys.path.append("../../utils")
import aou_utils


def UploadGS(local_path, gcp_path):
	"""
	Upload a local file to GCP

	Arguments
	---------
	local_path : str
	   Local path
	gcp_path : str
	   GCP path to upload to
	"""
	cmd = "gsutil cp {src} {dest}".format(src=local_path, dest=gcp_path)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))	

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the job", required=True, type=str)
	parser.add_argument("--tr", help="List of TRs to be extracted", required=True, type=str)
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
	
    # Upload TR bed file
	if args.tr.startswith("gs://"):
		tr_gcs = args.tr
	else:
		tr_gcs = gs_prefix + args.name + "/" + args.tr
		UploadGS(args.tr, tr_gcs)

	# Set up workflow JSON
	json_dict = {}
	json_dict["tr_extraction.vcfs"] =  [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.vcf.gz')]
	json_dict["tr_extraction.vcfs_index"] =  [gs_prefix + blob.name for blob in bucket.list_blobs(prefix=pfile) if blob.name.endswith('.vcf.gz.tbi')]
	json_dict["tr_extraction.str"] = args.tr
	json_dict["tr_extraction.out_prefix"] = args.name
	

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
	aou_utils.RunWorkflow("extract_str.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
