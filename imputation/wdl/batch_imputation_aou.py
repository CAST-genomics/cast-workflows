#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code to impute 10 samples at CBL region 

./batch_imputation_aou.py \
--name batch_test 
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.bref3 \
--mem 120 \
--batch-num 2
"""


import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import csv
from utils import MSG, ERROR

def GetFileBatches(sample_list, batch_num=None):
    sample_batch = []
    
    # Open txt file
    with open(sample_list, "r") as f:
        if batch_num is None:
            # Read all lines if batch_num is None
            sample_batch = f.readlines()
        else:
            # Read the first batch_num lines
            for _ in range(batch_num):
                line = f.readline().strip()  # Read and strip newline character
                if line:  # Check if the line is not empty
                    sample_batch.append(line)
                else:
                    break  # Exit the loop if there are no more lines
    
    return sample_batch

def RunWorkflow(json_file, json_options_file, wdl_dependencies_file, cromwell, dryrun=False):
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
		cmd = "cromshell submit ../wdl/batch_imputation.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	else:
		cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
	  			"cromwell-87.jar run imputation.wdl " + \
	  			"--inputs {} --options {}".format(json_file, json_options_file)
		
	if wdl_dependencies_file.strip() != "":
		cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))


def DownloadGS(filename):
	"""
	Download a GCP path locally

	Arguments
	---------
	filename : str
	   GCP path
	"""
	cmd = "gsutil cp {filename} .".format(filename=filename)
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

def ZipWDL(wdl_dependencies_file):
	"""
	Put all WDL dependencies into a zip file

	Arguments
	---------
	wdl_dependencies_fie : str
	    Zip file to put other wdls in
	"""
	files = ["imputation.wdl","processTR.wdl","processSNP.wdl","merge_TR_batch.wdl","merge_SNP_batch.wdl"]
	dirname = tempfile.mkdtemp()
	for f in files:
		shutil.copyfile("../wdl/%s"%f, dirname+"/"+f)
	shutil.make_archive(os.path.splitext(wdl_dependencies_file)[0], "zip", root_dir=dirname)


def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--vcf", help="Name of the genotype vcf file", required=True, type=str)
	parser.add_argument("--ref-panel", help="File id of ref genome", type=str)
	#parser.add_argument("--ref", help="File of ref genome", type=str,required=True)
	#parser.add_argument("--ref-index", help="File of ref index genome", type=str,required=False)
	parser.add_argument("--mem", help="Specify run memory ", type=int, required=False, default=50)
	parser.add_argument("--disk", help="Specify disk memory ", type=int, required=False, default=25)
	parser.add_argument("--window", help="Specify window size for imputation ", type=int, required=False, default=20)
	parser.add_argument("--samples", help="List of samples to process ", type=str, required=False, \
					 default="gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/tr_imputation/tr_imputation/sample/sample_manifest.txt")
	parser.add_argument("--region", help="Name of chrom position  chr:xxx-xxx", type=str,required=False)
	parser.add_argument("--beagle-region", help="Apply chrom for beagle", action="store_true",required=False)
	parser.add_argument("--subset-region", help="Subsetting region for vcf file", action="store_true",required=False)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)",type=int, required=False, default=None)
	parser.add_argument("--overlap", help="Specify overlap size for imputation ", type=int, required=False, default=2)
	parser.add_argument("--map", help="Specify genetic map for imputation ", type=str, required=True)					
	#parser.add_argument("--header-file", help="Add hipstr header",type=str, required=False, \
	#				 default="gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/tr_imputation/tr_imputation/header_annotation.txt")
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

	#Set up sample file list
	if args.samples.startswith("gs://"):
		DownloadGS(args.samples)
		sample_list = os.path.basename(args.samples)
	else: sample_list = args.samples

	# Upload vcf file
	if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)

	# set up batches of file
	sample_batch = GetFileBatches(sample_list,args.batch_num)

	# Upload subset sample file
	#if args.samples.startswith("gs://"):
	#	sample_file = args.samples
	#else:
	#			# Copying the exclude sample file
	#	sample_file = output_bucket + "/" + args.name + "/"
	#	UploadGS(args.samples, sample_file)


	if args.subset_region and args.region is None:
		ERROR("Must specify --region for --subset-region")



	# Set up workflow JSON
	json_dict = {}
	json_dict["batch_imputation.vcf"] = args.vcf
	json_dict["batch_imputation.vcf_index"]=args.vcf+".tbi"
	json_dict["batch_imputation.ref_panel"] = args.ref_panel
	json_dict["batch_imputation.out_prefix"] = args.name
	json_dict["batch_imputation.GOOGLE_PROJECT"] = project
	json_dict["batch_imputation.mem"] = args.mem
	json_dict["batch_imputation.disk"] = args.disk
	json_dict["batch_imputation.window_size"] = args.window
	json_dict["batch_imputation.samples"] = sample_batch 
	json_dict["batch_imputation.region"] = args.region
	json_dict["batch_imputation.subset_region"] = args.subset_region 
	json_dict["batch_imputation.beagle_region"] =args.beagle_region
	#json_dict["batch_imputation.header_file"] =args.header_file
	#json_dict["batch_imputation.ref"] =args.ref
	#json_dict["batch_imputation.ref_index"] =args.ref + ".tbi"
	json_dict["batch_imputation.overlap"] =args.overlap
	json_dict["batch_imputation.map"] =args.map




	


	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json optionsß
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Zip all the WDL depencies
	wdl_dependencies_file = args.name + "-wdl.zip"
	ZipWDL(wdl_dependencies_file)


	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file, wdl_dependencies_file,args.cromwell, dryrun=args.dryrun)


if __name__ == "__main__":
	main()
