#!/usr/bin/env python3
"""
./process_cram_list.py ukb_cram_files_long.txt > ukb_cram_and_index_files.txt

Output file of:
<file ID cram> <file ID index>
"""

import sys
import pandas as pd
# Problematic samples in 200k UKB release
# rmsamples = [
#   "1041104_23193_0_0",
#   "1523211_23193_0_0",
#   "2378461_23193_0_0",
#   "5172697_23193_0_0",
#   "5018133_23193_0_0",
#   "4166179_23193_0_0",
#   "4135688_23193_0_0",
#   "2084527_23193_0_0",
#   "4254386_23193_0_0",
#   "2379887_23193_0_0"
# ]
# For 500k samples, start with all samples
qced_sample = pd.read_csv("./QCed_sampleID_for_ALL.txt", sep="\t", header=None, names=["ID"])
# 500k cram file prfix 1000011_24048_0_0
keepsamples = qced_sample["ID"].apply(lambda x: f"{x}_24048_0_0").tolist()
#keepsamples = qced_sample["ID"].apply(lambda x: f"{x}_23372_0_0").tolist()

rm_samples = pd.read_csv("./samples_with_bad_cram.txt", sep="\t", header=None, names=["ID"])
rmsamples = rm_samples["ID"].tolist()

fdict = {} # sample -> {cram, idx}

with open(sys.argv[1], "r") as f:
    for line in f:
        if not line.startswith("closed"): continue
        items = line.strip().split()
        fname = items[5]
        fid = items[6].replace("(","").replace(")","")
        sample = fname.split(".")[0]
        if sample not in fdict:
            fdict[sample] = {}
        if fname.endswith(".crai"):
            fdict[sample]["idx"] = fid
        else:
            fdict[sample]["cram"] = fid

for sample in fdict.keys():
    if sample not in keepsamples: continue
    if sample in rmsamples: continue
    cramid = fdict[sample].get("cram", None)
    idxid = fdict[sample].get("idx", None)
    if cramid is None or idxid is None:
        sys.stderr.write("Couldn't find both files for %s"%sample)
    else:
        sys.stdout.write(" ".join([cramid, idxid])+"\n")
