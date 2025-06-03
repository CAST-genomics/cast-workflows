#!/usr/bin/env python3
"""
Script to rerun annotaTR AOU TR imputation to output new INFO field DSCOUNT which contains the non-zero dosage and their counts number in the pvar file.


chrom=1  
./rerun_annotaTR_aou.py \
--name chr${chrom}_annotated \
--chrom ${chrom} \

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
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--chrom", help="Which chromosome to process", type=str, required= True)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	args = parser.parse_args()
	

	# Set up workflow JSON
	json_dict = {}
	json_dict["rerun_annotator.vcf"] = os.environ.get("WORKSPACE_BUCKET") + "/imputation_merged_TR/chr%s_TR_merged.vcf.gz"%args.chrom
	json_dict["rerun_annotator.vcf_index"] = os.environ.get("WORKSPACE_BUCKET") + "/imputation_merged_TR/chr%s_TR_merged.vcf.gz.tbi"%args.chrom
	json_dict["rerun_annotator.ref_vcf"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.vcf.gz"%args.chrom
	json_dict["rerun_annotator.ref_index"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.vcf.gz.tbi"%args.chrom
	json_dict["rerun_annotator.out_prefix"] = args.name


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
	aou_utils.RunWorkflow("wdl/rerun_annotator.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
