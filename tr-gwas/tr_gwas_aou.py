#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 

./tr_gwas_aou.py 

"""

# need to change GetFileBatches
import argparse
import json
import os
import subprocess
import sys
import glob

sys.path.append("../utils")
import aou_utils

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--pgens", help="Path to pgen files", required=True, type=str)
	parser.add_argument("--psams", help="Path to psam files", required=True, type=str)
	parser.add_argument("--pvars", help="Path to pvar files", required=True, type=str)
	parser.add_argument("--phenotypes", help="Path to phenotype files",trequired=True, type=str)
	parser.add_argument("--cohorts", help="Path to cohort files",required=True, type=str)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	args = parser.parse_args()


	bucket = os.environ.get("WORKSPACE_BUCKET")
	# Set up workflow JSON
	json_dict = {}
	json_dict["tr_gwas.pgens"] = glob.glob(os.path.join(bucket,"tr_imputation/enstr-v3/results-250K/" + "*.pgen"))
	json_dict["tr_gwas.psams"] = glob.glob(os.path.join(bucket,"tr_imputation/enstr-v3/results-250K/" + "*.psam"))
	json_dict["tr_gwas.pvars"] = glob.glob(os.path.join(bucket,"tr_imputation/enstr-v3/results-250K/" + "*.pvar"))
	json_dict["tr_gwas.phenotypes"] = os.environ.get("WORKSPACE_BUCKET") + "/phenotypes/"
	json_dict["tr_gwas.cohorts"] = os.environ.get("WORKSPACE_BUCKET") + "/samples/"
	

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
	aou_utils.RunWorkflow("wdl/tr_gwas.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)



if __name__ == "__main__":
	main()
