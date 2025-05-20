#!/usr/bin/env python3

import json
import os

json_dict = json.loads(open("chr9concat.aou.json", "r").read())

# Replace last file with fixed one
for flist in json_dict["concatenate_batch_vcfs.batch_vcf_files"]:
	batchnum = flist[0].split("/")[-1].split("-")[0].replace("batch","")
	flist[-1] = os.environ["WORKSPACE_BUCKET"] + \
		"/subset_vcf/fix_chr9/" + "batch%s-chr9-130000000-138394717.vcf.gz"%batchnum

with open("chr9concat.fix.aou.json", "w") as f:
	json.dump(json_dict, f, indent=4)

