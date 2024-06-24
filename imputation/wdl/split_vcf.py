#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 
Example command:  python split_vcf.py --name extract_sample_aou --vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz
output file:  ${WORKSPACE_BUCKET}/cromwell-execution/split_vcf/d2b62462-a509-4a78-b418-b0dd3b04667e/call-get_sample_list/extract_sample_aou.txt

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
	cmd = "cromshell submit ../wdl/split_vcf.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	#cmd = "java -jar -Dconfig.file={} ".format("cromwell.conf") + \
	#			"cromwell-86.jar run split_vcf.wdl " + \
	#			"--inputs {} --options {}".format(json_file, json_options_file)
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

	# Upload vcf file
	if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = args.name + "/"
		#vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)

	# Set up workflow JSON
	json_dict = {}
	json_dict["split_vcf.vcf"] = args.vcf
	json_dict["split_vcf.vcf_index"]= args.vcf+".tbi"
	json_dict["split_vcf.out_prefix"] = args.name
	json_dict["split_vcf.GOOGLE_PROJECT"] = project
	json_dict["split_vcf.GCS_OAUTH_TOKEN"] = token

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
	RunWorkflow(json_file, json_options_file, dryrun=args.dryrun)


if __name__ == "__main__":
	main()
