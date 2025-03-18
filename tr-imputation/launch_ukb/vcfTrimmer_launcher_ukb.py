#!/usr/bin/env python3
"""
Script to launch UKB targeted TR analysis
"""

import argparse
import dxpy
import sys
import time
import json
def GetFileBatches(file_list, batch_size, batch_num=-1):
	"""
	Given a file containing a list of file IDs for the 
	vcfs, divide into batches of 
	size batch_size. 

	Arguments
	---------
	file_list : str
	   File containing a list of IDs for vcfs must have 
       one row per vcf, and the file IDs must be 
       separated by whitespace
	batch_size : int
	   Number of vcf files to include in each batch
	batch_num : int
	   Number of batches to create. If -1, make batches
	   out of everything. Primarily used for debugging
	   to make small test sets

	Returns
	-------
	vcf_batches : list of lists
	   Each internal list is a list of file IDs for pVCF files 
	   To play nicely with dxCompiler json, each item
	   in the internal list has format {"$dnanexus_link": pVCF_id}
	"""
	# Keep track of batches
	vcf_batches = []

	# VCFs in current batch
	curr_batch_vcfs = []
	with open(file_list, "r") as f:
		for line in f:
			# If we reached the batch size, add batch to 
			# the list and start a new one
			if len(curr_batch_vcfs) == batch_size:
				vcf_batches.append(curr_batch_vcfs)
				# Quit if we reached our desired number of batches
				# If -1 we just keep going until we run out of files
				if len(vcf_batches) == batch_num: break 
				curr_batch_vcfs = []
			# Add this vcf to the current batch
			vcf_id= line.strip()
			curr_batch_vcfs.append(dxpy.dxlink(vcf_id))
	# Add any leftovers
	if len(curr_batch_vcfs) > 0 and \
		(len(vcf_batches)<batch_num or batch_num==-1):
		vcf_batches.append(curr_batch_vcfs)
	return vcf_batches 

def saveJson(file_name, input_dict):
    with open(file_name, "w") as file:
        json.dump(input_dict, file, indent=4)

def RunWorkflow(json_dict, workflow_id, name, depends=[]):
	"""
	Run DNA Nexus workflow

	Arguments
	---------
	json_dict : dict
		Input arguments for the workflow
	workflow_id : str
		ID of the DNA Nexus worfklow
	name : str
		Used to determine where to store output

	Returns
	-------
	analysis : DXAnalysis object
	"""
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
	analysis = workflow.run(json_dict, \
		depends_on=[item.get_id() for item in depends], \
		folder="/yal084/tr_imputation/trimmed_vcfs/{name}".format(name=name))
	sys.stderr.write("Started analysis %s (%s)\n"%(analysis.get_id(), name))
    # sleep to delay launches between analyses as each analysis mounts
    # a bunch of files which can put pressure on the API call limits
    # imposed by DNANexus that they care about but don't cleanly enforce
	time.sleep(30)
	return analysis

# def UploadDNANexus(fname, name):
# 	"""
# 	Upload a file to DNA Nexus and return
# 	its file ID
#
# 	Arguments
# 	---------
# 	fname : str
# 	   Path to file to upload
# 	name : str
# 	   Name of subdirectory to store in
#
# 	Returns
# 	-------
# 	file : dxlink
# 	   {"$dnanexus_link": file_id}
# 	"""
# 	sys.stderr.write("Uploading BED file...\n")
# 	folder = "/TargetedSTR/results/{name}".format(name=name)
# 	dxfile = dxpy.upload_local_file(fname, folder=folder, parents=True)
# 	print(dxpy.dxlink(dxfile))
# 	return dxpy.dxlink(dxfile)

def GetJBOR(analysis, filename):
	return {
		"$dnanexus_link": {
			"analysis": analysis.get_id(),
			"field": filename
		}
	}

def main():
	parser = argparse.ArgumentParser(__doc__)
#	parser.add_argument("--tr-bed", help="BED file with TR regions to genotype", required=True, type=str)
	parser.add_argument("--chrom", help="chromosome of current pVCFs e.g. chr1", required=True, type=str)
	parser.add_argument("--batch-size", help="trimed vcf batch size", type=int, default=200)
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)", required=False, default=-1, type=int)
	parser.add_argument("--qc-thresholds", help="threshold for qc", type=str, required=False, default="INFO/AAScore>=0.5")
	parser.add_argument("--rm-tags", help="tags to remove", required=False, default="FORMAT/FT,FORMAT/AD,FORMAT/MD,FORMAT/DP,FORMAT/RA,FORMAT/PP,FORMAT/GQ", type=str)
	parser.add_argument("--threads-num", help="num of threads for bcftools", required=False, default=1, type=int)

	# parser.add_argument("--file-list", help="List of vcfs to process"
	# 	"Format of each line: vcf-file-id", type=str, required=False, default="ukb_pVCF_list.txt")
	# parser.add_argument("--genome-id", help="File id of ref genome", type=str, default="file-GGJ1z28JbVqbpqB93YbPqbzz")
	# parser.add_argument("--genome-idx-id", help="File id of ref genome index", type=str, default="file-GGJ94JQJv7BGFYq8BGp62xPV")
	parser.add_argument("--workflow-id", help="DNA Nexus workflow ID", required=False, default="workflow-GzGv6qQJX3JX371q98yJ61fG")
	# Options for multi-batches
	parser.add_argument("--concat-workflow-id", help="DNA Nexus workflow ID for merging", required=False, default="workflow-GzGJ6b8JX3JjYZVZ3ZY83P4z")
	parser.add_argument("--max-batches-per-workflow", help="Maximum number of batches to launch at once. -1 means all", required=False, default=10, type=int)
	parser.add_argument("--concurrent", help="Launch all batches at once", action="store_true")
	parser.add_argument("--bcftools-mem", help="Bcftools run memory, modify if run smaller sample size, default=16G", required=False, type=int, default=16)
	parser.add_argument("--bcftools-threads", help="Bcftools threads, more threads requires large memory", required=False, type=int, default=1)
	parser.add_argument("--concat-mem", help="Merge hipstr run memory, modify if run smaller sample size, default=4G", required=False, type=int, default=16)
	args = parser.parse_args()

	# Set up workflow JSON
	json_dict = {}
	# json_dict["stage-common.genome"] = {}
	# json_dict["stage-common.genome_index"] = {}
	# json_dict["stage-common.genome"] = dxpy.dxlink(args.genome_id)
	# json_dict["stage-common.genome_index"] = dxpy.dxlink(args.genome_idx_id)
	json_dict["stage-common.outprefix"] = args.chrom
	json_dict["stage-common.bcftools_mem"] = args.bcftools_mem
	json_dict["stage-common.qc_thresholds"] = args.qc_thresholds
	json_dict["stage-common.rm_tags"] = args.rm_tags
	json_dict["stage-common.threads_num"] = args.threads_num

	# json_dict["stage-common.merge_mem"] = args.merge_mem
	# json_dict["stage-common.infer_samps_from_file"] = True

	# Upload
	# if args.tr_bed.startswith("file-"):
	# 	json_dict["stage-common.tr_bed"] = args.tr_bed
	# else:
	# 	json_dict["stage-common.tr_bed"] = UploadDNANexus(args.tr_bed, args.chrom)

	# Set up batches of files
	sys.stderr.write("Setting up batches...\n")
	vcf_batches = GetFileBatches(f"{args.chrom}_pvcf_file_ids.txt", args.batch_size, args.batch_num)

	# Run batches
	# final_vcf = None
	# final_vcf_idx = None
	if args.max_batches_per_workflow == -1:
		# Run all at once
		json_dict["stage-common.pvcf_batches"] = vcf_batches
#		json_dict["stage-common.cram_index_batches"] = cram_idx_batches
		analysis = RunWorkflow(json_dict, args.workflow_id, args.chrom)
	else:
		# Run in chunks across multiple workflows
		depends = []
		batch_num = 0
		curr_idx = 0
		curr_vcf_batches = []
#		curr_idx_batches = []
		while curr_idx < len(vcf_batches):
			if len(curr_vcf_batches) == args.max_batches_per_workflow:
				batch_name = args.chrom + "-CHUNK%s"%batch_num
				batch_dict = json_dict.copy()
				batch_dict["stage-common.outprefix"] = batch_name
				batch_dict["stage-common.pvcf_batches"] = curr_vcf_batches
#				batch_dict["stage-common.cram_index_batches"] = curr_idx_batches
				if args.concurrent:
					use_dep = []
				else: use_dep = depends
				analysis = RunWorkflow(batch_dict, args.workflow_id, \
					args.chrom + "/" + args.chrom+"_%s"%batch_num, depends=use_dep)
				depends.append(analysis)
				saveJson(f"./logs/{args.chrom}_{batch_num}.json", batch_dict)
				batch_num += 1
				curr_vcf_batches = []
#				curr_idx_batches = []
			curr_vcf_batches.append(vcf_batches[curr_idx])
#			curr_idx_batches.append(cram_idx_batches[curr_idx])
			curr_idx += 1
		# Run last batch if any are left
		if len(curr_vcf_batches) > 0:
			batch_name = args.chrom + "-CHUNK%s"%batch_num
			batch_dict = json_dict.copy()
			batch_dict["stage-common.outprefix"] = batch_name
			batch_dict["stage-common.pvcf_batches"] = curr_vcf_batches
#			batch_dict["stage-common.cram_index_batches"] = curr_idx_batches
			if args.concurrent:
				use_dep = []
			else: use_dep = depends
			analysis = RunWorkflow(batch_dict, args.workflow_id, \
				args.chrom + "/" + args.chrom+"_%s"%batch_num, depends=use_dep)
			depends.append(analysis)
			saveJson(f"./logs/{args.chrom}_{batch_num}.json", batch_dict)

		# Run a final job to merge all the meta-batches
		merge_vcfs = []
		merge_vcfs_idx = []
		for analysis in depends:
			vcf = GetJBOR(analysis, "stage-outputs.final_vcf")
			vcf_idx = GetJBOR(analysis, "stage-outputs.final_vcf_indx")
			merge_vcfs.append(vcf)
			merge_vcfs_idx.append(vcf_idx)
		sys.stderr.write("Setting up merge from meta-batches...\n")
		merge_dict = {}
		merge_dict["stage-common.merged_prefix"] = args.chrom
		merge_dict["stage-common.vcf_files"] = merge_vcfs
		merge_dict["stage-common.vcf_indexs"] = merge_vcfs_idx
		merge_dict["stage-common.mem"] = args.concat_mem 
		analysis = RunWorkflow(merge_dict, args.concat_workflow_id, args.chrom, depends=depends)
		saveJson(f"./logs/{args.chrom}_concat.json", merge_dict)
if __name__ == "__main__":
	main()
