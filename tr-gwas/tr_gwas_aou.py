#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 

./tr_gwas_aou.py 
--name 

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
	parser.add_argument("--name", help="name of the run", required=True, type=str)
	#parser.add_argument("--pgens", help="Path to pgen files", required=True, type=str)
	#parser.add_argument("--psams", help="Path to psam files", required=True, type=str)
	#parser.add_argument("--pvars", help="Path to pvar files", required=True, type=str)
	#parser.add_argument("--phenotypes", help="Path to phenotype files",required=True, type=str)
	#parser.add_argument("--cohorts", help="Path to cohort files",required=True, type=str)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	args = parser.parse_args()

	
	# Set up workflow JSON
	json_dict = {}
	json_dict["tr_gwas.pgens"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/pgen_list.txt"
	json_dict["tr_gwas.psams"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/psam_list.txt"
	json_dict["tr_gwas.pvars"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/pvar_list.txt"
	json_dict["tr_gwas.phenotypes"] = os.environ.get("WORKSPACE_BUCKET") + "/phenotypes/phenotype_list.txt"
	json_dict["tr_gwas.cohorts"] = os.environ.get("WORKSPACE_BUCKET") + "/samples/sample_list.txt"
	
	

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
