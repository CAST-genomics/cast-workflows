import pandas as pd
import numpy as np
import os
import hail as hl
import subprocess
import sys

#default params
num_batches = 5 #recommend batches of 50k samples
batch_file_name = "wgs_ratios_batch"

if len(sys.argv) == 2:
    num_batches = sys.argv[1]
if len(sys.argv) == 3:
    batch_file_name = sys.argv[2]
#TODO add argument for output location

bucket = os.getenv('WORKSPACE_BUCKET')

hl.init(default_reference = "GRCh38")

mt_wgs_path = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
mt = hl.read_matrix_table(mt_wgs_path)
sampleIDs = mt.cols().to_pandas()
mt = mt.annotate_cols(sint = hl.int(mt.s))


qs = np.linspace(0,1, num_batches+1)
quantiles = sampleIDs.astype('int').quantile(qs)
quantiles.reset_index(inplace=True)

for i in range(num_batches):
    start = quantiles['s'][i]
    end = quantiles['s'][i+1]
    if i+1 == num_batches: #last round
        end += 1 #for last round, want to end to include last sample
    #filtered mt
    fmt = mt.filter_cols(((mt.sint >= start) & (mt.sint < end)))
    fmt = hl.sample_qc(fmt)
    fmt = fmt.filter_cols(fmt.sample_qc.call_rate >= .90, keep = True)

    ratios = fmt.select_cols(singleton = fmt.sample_qc.n_singleton,
            ins_del = fmt.sample_qc.r_insertion_deletion,
            ti_tv = fmt.sample_qc.r_ti_tv,
            het_hom = fmt.sample_qc.r_het_hom_var)
    tbl = ratios.cols()
    pddf = tbl.to_pandas()
    f_towrite = batch_file_name+str(i)+".csv"
    pddf.to_csv(f_towrite, sep=",", index=False)
    subprocess.run("gsutil cp "+f_towrite+" "+bucket+"/tmirmira/data/qc/"+f_towrite, shell=True, check=True)

