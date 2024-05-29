#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 
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
	#cmd = "cromshell submit ../wdl/beagle.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
				"cromwell-87.jar run create_reference.wdl " + \
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
	parser.add_argument("--snp-vcf", help="Name of the genotype vcf file including snps", required=True, type=str)
	parser.add_argument("--region", help="region in the SNP file to be extracted", required=True, type=str)
	parser.add_argument("--vntr-vcf", help="Name of the genotype vcf file including vntr(s)", required=True, type=str)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")

	args = parser.parse_args()

	
	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name
	output_path = os.path.join(bucket, "workflows", "cromwell-executions", "vntr_ref")

	# Upload vcf file
	'''if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)
        '''
	# Set up workflow JSON
	json_dict = {}
	json_dict["create_reference.snp_vcf"] = args.snp_vcf
	json_dict["create_reference.snp_vcf_index"]=args.snp_vcf+".tbi"
	json_dict["create_reference.vntr_vcf"] = args.vntr_vcf
	json_dict["create_reference.vntr_vcf_index"]=args.vntr_vcf+".tbi"
	json_dict["create_reference.region"] = args.region
	json_dict["create_reference.out_prefix"] = args.name
	json_dict["create_reference.GOOGLE_PROJECT"] = project

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json options
	json_options_dict = {"jes_gcs_root": output_path}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file, dryrun=args.dryrun)


if __name__ == "__main__":
	main()
