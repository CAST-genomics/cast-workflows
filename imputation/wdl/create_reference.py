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


def RunWorkflow(json_file, json_options_file, wdl_dependencies_file, workflow, cromwell=False, dryrun=False):
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
	if not cromwell:
	    cmd = "cromshell submit {workflow}.wdl {json} -op {options}".format(
                        workflow=workflow, json=json_file, options=json_options_file)
	else:
	    cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
				"cromwell-87.jar run {}.wdl ".format(workflow) + \
				"--inputs {} --options {}".format(json_file, json_options_file)
	if wdl_dependencies_file.strip() != "":
		cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
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
	files = ["create_reference.wdl"]
	dirname = tempfile.mkdtemp()
	for f in files:
		shutil.copyfile("../wdl/%s"%f, dirname+"/"+f)
	shutil.make_archive(os.path.splitext(wdl_dependencies_file)[0], "zip", root_dir=dirname)

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--snp-vcf", help="Name of the genotype vcf file including snps", required=True, type=str)
	parser.add_argument("--regions", help="A file with regions in the SNP file to be extracted",
                                required=True, type=str)
	parser.add_argument("--samples", help="A file with Sample ids to be used for creating the reference",
                                required=True, type=str)
	parser.add_argument("--vntr-vcf", help="Name of the genotype vcf file including vntr(s)", required=True, type=str)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--cromwell", help="Run in cromwell instead of the default cromshell.", action="store_true")
	parser.add_argument("--mem", help="Memory and disk in each task", type=int, required=False, default=20)
	parser.add_argument("--window", help="Beagle window size", type=int, required=False, default=20)
	parser.add_argument("--map", help="Path to the genetic map file", required=True, type=str)
	parser.add_argument("--chrom", help="Specify target chromosome ", type=str, required=False)

	args = parser.parse_args()

	
	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name
	output_path = os.path.join(bucket, "workflows", "cromwell-executions", "vntr_ref_merge")
	workflow="create_reference"

	# Set up workflow JSON
	json_dict = {}
	json_dict[workflow + ".snp_vcf"] = args.snp_vcf
	json_dict[workflow + ".snp_vcf_index"]=args.snp_vcf+".tbi"
	json_dict[workflow + ".vntr_vcf"] = args.vntr_vcf
	json_dict[workflow + ".vntr_vcf_index"]=args.vntr_vcf+".tbi"
	json_dict[workflow + ".regions"] = args.regions
	json_dict[workflow + ".samples"] = args.samples
	json_dict[workflow + ".bucket"] = bucket
	json_dict[workflow + ".mem"] = args.mem
	json_dict[workflow + ".window"] = args.window
	if args.chrom is not None:
	    json_dict[workflow + ".chrom"] = args.chrom
	    json_dict[workflow + ".map"] = os.path.join(args.map, "plink.{}_w_chr.GRCh38.map".format(args.chrom))
	else:
	    json_dict[workflow + ".map"] = args.map
	json_dict[workflow + ".out_prefix"] = args.name
	json_dict[workflow + ".GOOGLE_PROJECT"] = project

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json options
	json_options_dict = {"jes_gcs_root": output_path}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Zip all the WDL depencies
	if workflow == "create_reference_batch":
		wdl_dependencies_file = "create_ref_batch_depenency_wdl.zip"
		ZipWDL(wdl_dependencies_file)
	else:
		wdl_dependencies_file = ""

	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file,
                    wdl_dependencies_file=wdl_dependencies_file,
                    workflow=workflow,
                    dryrun=args.dryrun, cromwell=args.cromwell)


if __name__ == "__main__":
	main()
