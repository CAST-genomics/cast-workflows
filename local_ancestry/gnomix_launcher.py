#!/usr/bin/env python3
"""
Script to launch AOU Gnomix for local ancestry inference

Dryrun:
chrom=11; ./gnomix_launcher.py \
  --name test-chr${chrom} \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom} \
  --max-batches 2 \
  --dryrun
"""

import argparse
import glob
import json
import os
import subprocess
import sys

sys.path.append("../utils")
import aou_utils

GNOMIXMODEL = "gs://artifacts.ucsd-medicine-cast.appspot.com/resources/pretrained_gnomix_models.tar.gz"

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--vcfdir", help="GCP bucket with batch VCF files", \
		type=str, required=True)
	parser.add_argument("--name", help="Name of the run", type=str, required=True)
	parser.add_argument("--chrom", help="Which chromosome to process", type=str, required=True)
	parser.add_argument("--gnomix-model", help="GCP path to pretrained model", type=str, \
		default=GNOMIXMODEL)
	# For debugging
	parser.add_argument("--dryrun", help="Set up job but don't launch", \
		action="store_true")
	parser.add_argument("--max-batches", help="Only run this many batches", \
		type=int, default=-1)
	args = parser.parse_args()

	batch_vcf_files = aou_utils.GetBatchVCFFiles(args.vcfdir, args.max_batches)

	# Set up workflow JSON
	json_dict = {}
	json_dict["local_ancestry.out_prefix"] = args.name
	json_dict["local_ancestry.batch_vcf_files"] = batch_vcf_files
	json_dict["local_ancestry.chrom"] = args.chrom
	json_dict["local_ancestry.gnomix_model"] = args.gnomix_model
	json_dict["local_ancestry.chainfile"] = os.environ.get("WORKSPACE_BUCKET") + "/gnomix/resources/hg38ToHg19.over.chain.gz"
	if args.chrom in ["8", "12", "14", "17"]:
		# Some chroms had duplicate entries in the 1000G phased VCF :(
		json_dict["local_ancestry.refpanel"] = os.environ.get("WORKSPACE_BUCKET") + "/gnomix/resources/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-dedup.vcf.gz"%args.chrom
		json_dict["local_ancestry.refpanel_index"] = os.environ.get("WORKSPACE_BUCKET") + "/gnomix/resources/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-dedup.vcf.gz.tbi"%args.chrom
	else:
		json_dict["local_ancestry.refpanel"] = os.environ.get("WORKSPACE_BUCKET") + "/gnomix/resources/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"%args.chrom
		json_dict["local_ancestry.refpanel_index"] = os.environ.get("WORKSPACE_BUCKET") + "/gnomix/resources/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"%args.chrom

	# Convert to JSON and save to a file
	json_file = args.name + ".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Zip all the WDL depencies
	wdl_dependencies_file = args.name + "-wdl.zip"
	aou_utils.ZipWDL(wdl_dependencies_file, ["wdl/gnomix.wdl"])

	# Run workflow on AoU using cromshell
	aou_utils.RunWorkflow("wdl/local_ancestry.wdl", json_file, json_options_file, \
		wdl_dependencies_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()