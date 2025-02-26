#!/usr/bin/env python3
"""
Script to launch AOU phewas

example code to run phewas on the TMCO1 locus

python phewas_aou.py --locus "TMCO1" \
                     --chrom "chr1" \
                     --start 165761972 \
                     --n-threads 10 \
                     --tr-vcf $tr_vcf
"""


import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from math import ceil
import pandas as pd


def RunWorkflow(wdl_file, json_file, json_options_file, wdl_dependencies_file, cromwell, dryrun=False):
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
		cmd = "cromshell submit {wdl} {json} -op {options}".format(
                        wdl=wdl_file,
                        json=json_file, options=json_options_file)
		if wdl_dependencies_file.strip() != "":
	  		cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
	else:
		cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
	  			"jar_files/cromwell-87.jar run {} ".format(wdl_file) + \
	  			"--inputs {} --options {}".format(json_file, json_options_file)
	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))

def ZipWDL(wdl_dependencies_file):
	"""
	Put all WDL dependencies into a zip file

	Arguments
	---------
	wdl_dependencies_file : str
	    Zip file to put other wdls in
	"""
	files = ["imputation.wdl"]
	dirname = tempfile.mkdtemp()
	for f in files:
		shutil.copyfile("../wdl/%s"%f, dirname+"/"+f)
	shutil.make_archive(os.path.splitext(wdl_dependencies_file)[0], "zip", root_dir=dirname)

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

def run_single_batch(args, samples_files):
	wdl_workflow = "phewas"
	wdl_file = "./{}.wdl".format(wdl_workflow)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name
	output_path = os.path.join(bucket, "workflows", "cromwell-executions", "sara_phewas")


	# Set up workflow JSON
	json_dict = {}
	json_dict[wdl_workflow + ".chrom"] = args.chrom


	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json options
	json_options_dict = {"jes_gcs_root": output_path}
	#json_options_dict = {"call-caching.enabled": "false"}
	#json_options_dict = {"callCaching.allowResultReuse": "false"}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	RunWorkflow(wdl_file,
                    json_file,
                    json_options_file,
                    wdl_dependencies_file,
                    args.cromwell,
                    dryrun=args.dryrun)


def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--tr-vcf", help="Path to the genotyped or imputed tr-vcf file.", required=True, type=str)
	parser.add_argument("--chrom", help="The chromosome where the locus is located.", type=str, required=True)
	parser.add_argument("--start", help="Start position for the TR locus.", type=int, required=True)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--cromwell", help="Run using cormwell as opposed to the default cromshell",
                            action="store_true", default=False)


	args = parser.parse_args()




if __name__ == "__main__":
	main()
