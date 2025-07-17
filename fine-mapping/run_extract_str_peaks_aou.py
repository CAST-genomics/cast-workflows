#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 
Example:
./run_extract_str_peaks_aou.py \
--name hdl_chrom8_str_peaks \
--chrom 8 \
--region test_hdl_chrom8_region.txt 

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
	parser.add_argument("--chrom", help="Name of the chromsome file", required=True, type=str)
	parser.add_argument("--name", help="Name of the output file", required=True, type=str)
	parser.add_argument("--region", help="Name of the region file", type=str, required=True)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	args = parser.parse_args()


	# Set up workflow JSON
	json_dict = {}
	json_dict["extract_str_peak_gt.pgen"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/results-250K/chr%s_annotated.pgen"%args.chrom
	json_dict["extract_str_peak_gt.psam"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/results-250K/chr%s_annotated.psam"%args.chrom
	json_dict["extract_str_peak_gt.pvar"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/results-250K/chr%s_annotated.pvar"%args.chrom
	json_dict["extract_str_peak_gt.region"] = args.region
	json_dict["extract_str_peak_gt.out_prefix"] = args.name


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
	aou_utils.RunWorkflow("wdl/extract_str_peak_gt.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
