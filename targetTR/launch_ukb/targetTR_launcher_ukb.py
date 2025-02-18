#!/usr/bin/env python3
"""
Script to launch UKB targeted TR analysis
"""

import argparse
import dxpy
import sys
import time

def GetFileBatches(file_list, batch_size, batch_num=-1):
	"""
	Given a file containing a list of file IDs for the 
	crams and their indices, divide into batches of 
	size batch_size. 

	Arguments
	---------
	file_list : str
	   File containing a list of IDs for crams and indices
	   File must have one row per CRAM/index, and the file
	   IDs must be separated by whitespace
	batch_size : int
	   Number of cram files to include in each batch
	batch_num : int
	   Number of batches to create. If -1, make batches
	   out of everything. Primarily used for debugging
	   to make small test sets

	Returns
	-------
	cram_batches : list of lists
	   Each internal list is a list of file IDs for crams
	   To play nicely with dxCompiler json, each item
	   in the internal list has format {"$dnanexus_link": cram_id}
	cram_idx_batches : list of lists
	   Each internal list is a list of file IDs for cram indices
	   To play nicely with dxCompiler json, each item
	   in the internal list has format {"$dnanexus_link": cram_id}
	"""
	# Keep track of batches
	cram_batches = []
	cram_idx_batches = []

	# Crams/indices in current batch
	curr_batch_crams = []
	curr_batch_indices = []
	with open(file_list, "r") as f:
		for line in f:
			# If we reached the batch size, add batch to 
			# the list and start a new one
			if len(curr_batch_crams) == batch_size:
				cram_batches.append(curr_batch_crams)
				cram_idx_batches.append(curr_batch_indices)
				# Quit if we reached our desired number of batches
				# If -1 we just keep going until we run out of files
				if len(cram_batches) == batch_num: break 
				curr_batch_crams = []
				curr_batch_indices = []
			# Add this cram to the current batch
			cram_id, idx_id = line.strip().split()
			curr_batch_crams.append(dxpy.dxlink(cram_id))
			curr_batch_indices.append(dxpy.dxlink(idx_id))
	# Add any leftovers
	if len(curr_batch_crams) > 0 and \
		(len(cram_batches)<batch_num or batch_num==-1):
		cram_batches.append(curr_batch_crams)
		cram_idx_batches.append(curr_batch_indices)
	assert(len(cram_batches) == len(cram_idx_batches))
	return cram_batches, cram_idx_batches

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
		folder="/TargetedSTR/results/{name}".format(name=name))
	sys.stderr.write("Started analysis %s (%s)\n"%(analysis.get_id(), name))
    # sleep to delay launches between analyses as each analysis mounts
    # a bunch of files which can put pressure on the API call limits
    # imposed by DNANexus that they care about but don't cleanly enforce
	time.sleep(30)
	return analysis

def UploadDNANexus(fname, name):
	"""
	Upload a file to DNA Nexus and return
	its file ID

	Arguments
	---------
	fname : str
	   Path to file to upload
	name : str
	   Name of subdirectory to store in

	Returns
	-------
	file : dxlink
	   {"$dnanexus_link": file_id}
	"""
	sys.stderr.write("Uploading BED file...\n")
	folder = "/TargetedSTR/results/{name}".format(name=name)
	dxfile = dxpy.upload_local_file(fname, folder=folder, parents=True)
	print(dxpy.dxlink(dxfile))
	return dxpy.dxlink(dxfile)

def GetJBOR(analysis, filename):
	return {
		"$dnanexus_link": {
			"analysis": analysis.get_id(),
			"field": filename
		}
	}

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--tr-bed", help="BED file with TR regions to genotype", required=True, type=str)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--batch-size", help="HipSTR batch size", type=int, default=300)
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)", required=False, default=-1, type=int)
	parser.add_argument("--file-list", help="List of crams and indices to process"
		"Format of each line: cram-file-id cram-index-id", type=str, required=False, default="ukb_cram_and_index_files.txt")
	parser.add_argument("--genome-id", help="File id of ref genome", type=str, default="file-GGJ1z28JbVqbpqB93YbPqbzz")
	parser.add_argument("--genome-idx-id", help="File id of ref genome index", type=str, default="file-GGJ94JQJv7BGFYq8BGp62xPV")
	parser.add_argument("--workflow-id", help="DNA Nexus workflow ID", required=False, default="workflow-GpkxfqjJv7B4pGXPp5K9V6q2")
	# Options for multi-batches
	parser.add_argument("--merge-workflow-id", help="DNA Nexus workflow ID for merging", required=False, default="workflow-Gpkxg18Jv7B23bPjy26GF642")
	parser.add_argument("--max-batches-per-workflow", help="Maximum number of batches to launch at once. -1 means all", required=False, default=10, type=int)
	parser.add_argument("--concurrent", help="Launch all batches at once", action="store_true")
	parser.add_argument("--hipstr-mem", help="Hipstr run memory, modify if run smaller sample size, default=16G", required=False, type=int, default=16)
	parser.add_argument("--merge-mem", help="Merge hipstr run memory, modify if run smaller sample size, default=4G", required=False, type=int, default=4)
	args = parser.parse_args()

	# Set up workflow JSON
	json_dict = {}
	json_dict["stage-common.genome"] = {}
	json_dict["stage-common.genome_index"] = {}
	json_dict["stage-common.genome"] = dxpy.dxlink(args.genome_id)
	json_dict["stage-common.genome_index"] = dxpy.dxlink(args.genome_idx_id)
	json_dict["stage-common.outprefix"] = args.name
	json_dict["stage-common.hipstr_mem"] = args.hipstr_mem
	json_dict["stage-common.merge_mem"] = args.merge_mem
	json_dict["stage-common.infer_samps_from_file"] = True

	# Upload
	if args.tr_bed.startswith("file-"):
		json_dict["stage-common.tr_bed"] = args.tr_bed
	else:
		json_dict["stage-common.tr_bed"] = UploadDNANexus(args.tr_bed, args.name)

	# Set up batches of files
	sys.stderr.write("Setting up batches...\n")
	cram_batches, cram_idx_batches = GetFileBatches(args.file_list, args.batch_size, args.batch_num)

	# Run batches
	final_vcf = None
	final_vcf_idx = None
	if args.max_batches_per_workflow == -1:
		# Run all at once
		json_dict["stage-common.cram_file_batches"] = cram_batches
		json_dict["stage-common.cram_index_batches"] = cram_idx_batches
		analysis = RunWorkflow(json_dict, args.workflow_id, args.name)
	else:
		# Run in chunks across multiple workflows
		depends = []
		batch_num = 0
		curr_idx = 0
		curr_cram_batches = []
		curr_idx_batches = []
		while curr_idx < len(cram_batches):
			if len(curr_cram_batches) == args.max_batches_per_workflow:
				batch_name = args.name + "-CHUNK%s"%batch_num
				batch_dict = json_dict.copy()
				batch_dict["stage-common.outprefix"] = batch_name
				batch_dict["stage-common.cram_file_batches"] = curr_cram_batches
				batch_dict["stage-common.cram_index_batches"] = curr_idx_batches
				if args.concurrent:
					use_dep = []
				else: use_dep = depends
				analysis = RunWorkflow(batch_dict, args.workflow_id, \
					args.name + "/" + args.name+"_%s"%batch_num, depends=use_dep)
				depends.append(analysis)
				batch_num += 1
				curr_cram_batches = []
				curr_idx_batches = []
			curr_cram_batches.append(cram_batches[curr_idx])
			curr_idx_batches.append(cram_idx_batches[curr_idx])
			curr_idx += 1
		# Run last batch if any are left
		if len(curr_cram_batches) > 0:
			batch_name = args.name + "-CHUNK%s"%batch_num
			batch_dict = json_dict.copy()
			batch_dict["stage-common.outprefix"] = batch_name
			batch_dict["stage-common.cram_file_batches"] = curr_cram_batches
			batch_dict["stage-common.cram_index_batches"] = curr_idx_batches
			if args.concurrent:
				use_dep = []
			else: use_dep = depends
			analysis = RunWorkflow(batch_dict, args.workflow_id, \
				args.name + "/" + args.name+"_%s"%batch_num, depends=use_dep)
			depends.append(analysis)
		# Run a final job to merge all the meta-batches
		merge_vcfs = []
		merge_vcfs_idx = []
		for analysis in depends:
			vcf = GetJBOR(analysis, "stage-outputs.finalvcf")
			vcf_idx = GetJBOR(analysis, "stage-outputs.finalvcf_index")
			merge_vcfs.append(vcf)
			merge_vcfs_idx.append(vcf_idx)
		sys.stderr.write("Setting up merge from meta-batches...\n")
		merge_dict = {}
		merge_dict["stage-common.out_name"] = args.name
		merge_dict["stage-common.vcf_files"] = merge_vcfs
		merge_dict["stage-common.vcf_idxs"] = merge_vcfs_idx
		analysis = RunWorkflow(merge_dict, args.merge_workflow_id, args.name, depends=depends)

if __name__ == "__main__":
	main()
