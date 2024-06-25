#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code to impute 100 samples 
./run_imputation.py \
--name test_imputation 
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.bref3 \
--samples-file $WORKSPACE_BUCKET/tr_imputation/tr_imputation/subset_samples/aou_subset_samples_100.txt \
--region chr21:5101889-6101889
--subset_region
--mem 30
"""





import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from utils import MSG, ERROR


def RunWorkflow(json_file, json_options_file, cromwell, dryrun=False):
	"""
	Run workflow on AoU

	Arguments
	---------
	json_file : str
		JSON file path with input arguments
	json_options_file : str
		JSON with additional options for cromshell

	dryrun : bool
		Just print the command, don't actually run cromshell
	"""
	if cromwell is False:
		cmd = "cromshell submit ../wdl/imputation.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	else:
		cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
	  			"cromwell-87.jar run imputation.wdl " + \
	  			"--inputs {} --options {}".format(json_file, json_options_file)
	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))
	
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
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--vcf", help="Name of the genotype vcf file", required=True, type=str)
	parser.add_argument("--ref-panel", help="File id of ref genome", type=str)
	parser.add_argument("--mem", help="Specify run memory ", type=int, required=False, default=32)
	parser.add_argument("--window", help="Specify window size for imputation ", type=int, required=False, default=20)
	parser.add_argument("--sample-file", help="Name of sub_samples file ", type=str, required=True)
	parser.add_argument("--region", help="Name of chrom position  chr:xxx-xxx", type=str,required=False)
	parser.add_argument("--beagle-region", help="Apply chrom for beagle", action="store_true",required=False)
	parser.add_argument("--subset-region", help="Subsetting region for vcf file", action="store_true",required=False)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--cromwell", help="Run using cormwell as opposed to the default cromshell",
                            action="store_true", default=False)


	args = parser.parse_args()

	
	# Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name

	# Upload vcf file
	if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)

	# Upload subset sample file
	if args.sample_file.startswith("gs://"):
		sample_file = args.sample_file
	else:
				# Copying the exclude sample file
		sample_file = output_bucket + "/" + args.name + "/"
		UploadGS(args.sample_file, sample_file)


	if args.subset_region and args.region is None:
		ERROR("Must specify --region for --subset-region")



	# Set up workflow JSON
	json_dict = {}
	json_dict["run_imputation.vcf"] = args.vcf
	json_dict["run_imputation.vcf_index"]=args.vcf+".tbi"
	json_dict["run_imputation.ref_panel"] = args.ref_panel
	json_dict["run_imputation.out_prefix"] = args.name
	json_dict["run_imputation.GOOGLE_PROJECT"] = project
	json_dict["run_imputation.mem"] = args.mem
	json_dict["run_imputation.window_size"] = args.window
	json_dict["run_imputation.sample_file"] = args.sample_file 
	json_dict["run_imputation.region"] = args.region
	json_dict["run_imputation.subset_region"] = args.subset_region 
	


	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json optionsß
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file, args.cromwell, dryrun=args.dryrun)


if __name__ == "__main__":
	main()