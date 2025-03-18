#!/usr/bin/env python3
"""
Script to launch AOU TR imputation
 
chrom=21
./batch_imputation_aou.py \
--name chr${chrom}_batch_test \
--chrom ${chrom} \
--vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
--batch-num 2
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
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--vcfdir", help="GCP bucket with batch VCF files", type=str, required=True)
	parser.add_argument("--chrom", help="Which chromosome to process", type=str, required=True)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)",type=int, required=False, default=-1)
	args = parser.parse_args()

	batch_vcf_files = aou_utils.GetBatchVCFFiles(args.vcfdir, args.batch_num)

	# Set up workflow JSON
	json_dict = {}
	json_dict["batch_imputation.ref_vcf"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.vcf.gz"%args.chrom
	json_dict["batch_imputation.ref_index"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.vcf.gz.tbi"%args.chrom
	json_dict["batch_imputation.ref_panel"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/enstr-v3/ensembletr_refpanel_v3_chr%s.bref3"%args.chrom
	json_dict["batch_imputation.map"] = os.environ.get("WORKSPACE_BUCKET") + "/tr_imputation/tr_imputation/genetic_map/beagle_chr%s_b38.map"%args.chrom
	json_dict["batch_imputation.out_prefix"] = args.name
	json_dict["batch_imputation.batch_vcf_files"] = batch_vcf_files 

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
	aou_utils.RunWorkflow("wdl/batch_imputation.wdl", json_file, \
		json_options_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()
