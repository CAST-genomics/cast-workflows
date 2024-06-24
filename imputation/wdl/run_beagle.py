#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code:
./imputation_aou.py \
--name test_imputation 
--vcf $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/chr11_final_SNP_merged_additional_TRs.vcf.gz \
--mem 60
"""


import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile


def RunWorkflow(json_file, json_options_file, dryrun=False):
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
	cmd = "cromshell submit ../wdl/beagle.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	#cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
	#  			"cromwell-86.jar run beagle.wdl " + \
	#  			"--inputs {} --options {}".format(json_file, json_options_file)
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
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")


	args = parser.parse_args()

	
	# Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name
	output_path = os.path.join(bucket, "workflows", "cromwell-executions", "beagle_3")

	# Upload vcf file
	if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)

	
	# Set up workflow JSON
	json_dict = {}
	json_dict["beagle.vcf"] = args.vcf
	json_dict["beagle.vcf_index"]=args.vcf+".tbi"
	json_dict["beagle.ref_panel"] = args.ref_panel
	json_dict["beagle.ref_panel_index"] = args.ref_panel+".tbi"
	json_dict["beagle.out_prefix"] = args.name
	json_dict["beagle.GOOGLE_PROJECT"] = project
	json_dict["beagle.GCS_OAUTH_TOKEN"] = token
	json_dict["beagle.mem"] = args.mem
	json_dict["beagle.window_size"] = args.window

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json options
	json_options_dict = {"jes_gcs_root": output_path}
	json_options_file = args.name+".options.aou.json"
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file, dryrun=args.dryrun)


if __name__ == "__main__":
	main()
