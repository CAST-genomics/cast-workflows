#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 
Example:
./run_extract_str_peaks_aou.py \
--phenotye hdl \
--chrom 8 \
--tr-list test_hdl_chrom8_region

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
	parser.add_argument("--chrom", help="Name of the chromsome file", required=True, type=str)
	parser.add_argument("--tr-list", help="Name of the region file", type=str, required=True)
	parser.add_argument("--phenotype", help="Name of the phenotyoe", type=str, required=True)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")

	args = parser.parse_args()


	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)
	project = os.getenv("GOOGLE_PROJECT")
	bucket= os.getenv("WORKSPACE_BUCKET")
	

	# Upload region file
	if args.tr_list.startswith("gs://"):
		tr_list_gcs = args.tr_list
	else:
		tr_list_gcs = bucket + "/" + args.phenotype + "/" + args.tr_list
		UploadGS(args.tr_list,tr_list_gcs)



	# Set up workflow JSON
	json_dict = {}
	json_dict["extract_str_peak_gt.pgen"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/results-250K/chr%s_annotated.pgen"%args.chrom
	json_dict["extract_str_peak_gt.psam"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/results-250K/chr%s_annotated.psam"%args.chrom
	json_dict["extract_str_peak_gt.pvar"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/results-250K/chr%s_annotated.pvar"%args.chrom
	json_dict["extract_str_peak_gt.tr_list"] = tr_list_gcs 
	json_dict["extract_str_peak_gt.pheno"] = args.phenotype
	json_dict["extract_str_peak_gt.GOOGLE_PROJECT"] = project
	json_dict["extract_str_peak_gt.GCS_OAUTH_TOKEN"] = token


	# Convert to json and save as a file
	json_file = args.phenotype + args.chrom + "tr_finemapped.aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
	json_options_dict = {}
	json_options_file = args.phenotype + args.chrom + "tr_finemapped.options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	aou_utils.RunWorkflow("wdl/extract_str_peak_gt.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
