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

	batch_vcf_files = aou_utils.GetBatchVCFFiles(args.vcfdir, args.batch_num)

	# Set up workflow JSON
	json_dict = {}
	json_dict["convert_vcf_to_pgen.vcf"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.vcf.gz"%args.chrom
	json_dict["convert_vcf_to_pgen.vcf_index"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.vcf.gz.tbi"%args.chrom
	json_dict["convert_vcf_to_pgen.out_prefix"] = args.name


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
