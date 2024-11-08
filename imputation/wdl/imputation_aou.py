#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code to impute 10 samples at CBL region 
./imputation_aou.py \
--name test_imputation 
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr11.vcf.bgz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/chr11_final_SNP_merged_additional_TRs.vcf.gz \
--samples-file $WORKSPACE_BUCKET/tr_imputation/tr_imputation/aou_subset_samples_10.txt \
--regions-file $WORKSPACE_BUCKET/tr_imputation/tr_imputation/CBL_region.bed \
--mem 70
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
                        # There is no conf flag in cromshell
                        #conf="/home/jupyter/cromwell.conf")
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
	#wdl_workflow = "imputation"
	wdl_workflow = "imputation_batch"
	wdl_file = "./{}.wdl".format(wdl_workflow)

	# Get token. This is unnecessary as token in updated in the WDL workflow.
	#token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
	#	capture_output=True, check=True, encoding='utf-8')
	#token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name
	output_path = os.path.join(bucket, "workflows", "cromwell-executions", "sara_beagle_batch")

	if False:
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
            if args.samples_file.startswith("gs://"):
                    samples_file = args.samples_file
            else:
                                    # Copying the exclude sample file
                    samples_file = output_bucket + "/" + args.name + "/"
                    UploadGS(args.samples_file, samples_file)



            # Upload subset region file
            if args.regions_file.startswith("gs://"):
                    regions_file = args.regions_file
            else:
                                    # Copying the exclude sample file
                    regions_file = output_bucket + "/" + args.name + "/"
                    UploadGS(args.regions_file, regions_file)



	# Set up workflow JSON
	json_dict = {}
	json_dict[wdl_workflow + ".vcf"] = args.vcf
	json_dict[wdl_workflow + ".vcf_index"]=args.vcf+".tbi"
	json_dict[wdl_workflow + ".ref_panel_bref"] = args.ref_panel.replace(".vcf.gz", ".bref3")
	json_dict[wdl_workflow + ".ref_panel"] = args.ref_panel
	json_dict[wdl_workflow + ".ref_panel_index"] = args.ref_panel+".tbi"
	json_dict[wdl_workflow + ".genetic_map"] = args.map
	json_dict[wdl_workflow + ".out_prefix"] = args.name
	json_dict[wdl_workflow + ".GOOGLE_PROJECT"] = project
	json_dict[wdl_workflow + ".mem"] = args.mem
	json_dict[wdl_workflow + ".bucket"] = bucket
	json_dict[wdl_workflow + ".chrom"] = args.chrom
	json_dict[wdl_workflow + ".window_size"] = args.window
	json_dict[wdl_workflow + ".overlap"] = args.overlap
	json_dict[wdl_workflow + ".samples_files"] = samples_files
	json_dict[wdl_workflow + ".regions_file"] = args.regions_file


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


	# Zip all the WDL depencies
	wdl_dependencies_file = args.name + "-wdl.zip"
	ZipWDL(wdl_dependencies_file)


	# Run workflow on AoU using cromwell
	RunWorkflow(wdl_file,
                    json_file,
                    json_options_file,
                    wdl_dependencies_file,
                    args.cromwell,
                    dryrun=args.dryrun)

def get_all_sample_ids():
    # all samples are in "data/sample_ids.txt"
    # all samples passing the sample tests are in
    # data/passing_samples_v7.1.csv
    sample_df = pd.read_csv("data/passing_samples_v7.1.csv")
    return sample_df

def get_samples_file(samples_df, batch_idx, id_column_name="person_id"):
    # For batches that previously failed, change the input sample name,
    # So that the subset_vcf has to run again and do not copy the cached files
    # From a failed attempt.
    failed_batches_to_be_repeated = []
    if batch_idx in failed_batches_to_be_repeated:
        filename = "batch_{}_repeat_size_1000.txt".format(batch_idx)
    else:
        filename = "batch_{}_size_1000.txt".format(batch_idx)
    filepath = os.path.join("samples_files", filename)
    if not os.path.exists("samples_files"):
        os.mkdir("samples_files")
    with open(filepath, "w+") as samples_file:
        text = ""
        for idx, row in samples_df.iterrows():
            text += str(row[id_column_name])+"\n"
        samples_file.write(text)

    # Upload the samples batch file to the workspace
    bucket = os.getenv("WORKSPACE_BUCKET")
    gs_path = os.path.join(bucket, "saraj", "samples_files_batches", filename)
    # Do not upload the files, the files are there already.
    # Uncommend if batches are changed.
    #UploadGS(filepath, gs_path)
    return gs_path

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--vcf", help="Name of the genotype vcf file", required=True, type=str)
	parser.add_argument("--ref-panel", help="File id of ref genome", type=str)
	parser.add_argument("--mem", help="Specify run memory ", type=int, required=False, default=32)
	parser.add_argument("--map", help="Path to the genetic map file", required=True, type=str)
	parser.add_argument("--chrom", help="Specify target chromosome ", type=str, required=True)
	parser.add_argument("--window", help="Specify window size for imputation ", type=int, required=False, default=5)
	parser.add_argument("--overlap", help="Specify overlap size for imputation ", type=int, required=False, default=2)
	parser.add_argument("--samples-file", help="Name of sub_samples file ", type=str, required=True)
	parser.add_argument("--regions-file", help="Name of sub_region file ", type=str,required=False)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--cromwell", help="Run using cormwell as opposed to the default cromshell",
                            action="store_true", default=False)
	parser.add_argument("--batch-size", help="Number of samples in each batch ", type=int, required=False, default=1)


	args = parser.parse_args()

        # Get batch size and number
	sample_ids_df = get_all_sample_ids()
        # For small tests, grab a subset
	#sample_ids_df = sample_ids_df.iloc[range(args.batch_size*10)]

	num_samples = len(sample_ids_df)
	num_batches = ceil(num_samples / args.batch_size)

	samples_files = []
	failed_batches_to_avoid = []

	for batch_idx in range(num_batches):
            if batch_idx in failed_batches_to_avoid:
                continue
            # Set samples_files for each batch
            start = batch_idx * args.batch_size
            end = min(len(sample_ids_df), (batch_idx + 1) * args.batch_size)
            partial_samples_df = sample_ids_df.iloc[range(start, end)]
            samples_files.append(get_samples_file(partial_samples_df, batch_idx))

	run_single_batch(args, samples_files)


if __name__ == "__main__":
	main()
