#!/usr/bin/env python3

"""
Run GWAS on All of Us

Example:
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
"""

import argparse
from gwas_plotter import PlotManhattan, PlotQQ
import math
import os
import pandas as pd
import re
import sys
from utils import MSG, ERROR
import numpy as np
import scipy.stats as stats
import warnings
import hail as hl

GWAS_METHODS = ["hail"]
ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.1.csv")
DEFAULTSAMPLECSV = "passing_samples_v7.1.csv"

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def CheckRegion(region):
    if region is None: return True
    return re.match(r"\w+:\d+-\d+", region) is not None

def GetOutPath(phenotype, region, samplefile):
    method = 'segmented_hail'
    cohort = "ALL" #if samplefile == DEFAULTSAMPLECSV
    if samplefile != DEFAULTSAMPLECSV:
        cohort = samplefile[:-4] #chop off .csv at the end
    outprefix = "%s_%s_%s"%(phenotype, method, cohort)
    if region is not None:
        outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
    return outprefix + ".gwas"

def GetFloatFromPC(x):
    x = x.replace("[","").replace("]","")
    return float(x)

def LoadAncestry(ancestry_pred_path):
    if ancestry_pred_path.startswith("gs://"):
        if not os.path.isfile("ancestry_preds.tsv"):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(ancestry_pred_path))
        ancestry_pred_path = "ancestry_preds.tsv"
    ancestry = pd.read_csv(ancestry_pred_path, sep="\t")
    ancestry.rename({"research_id": "person_id"}, axis=1, inplace=True)
    num_pcs = len(ancestry["pca_features"].values[0].split(","))
    pcols = ["PC_%s"%i for i in range(num_pcs)]
    ancestry[pcols] = ancestry["pca_features"].str.split(",", expand=True)
    for p in pcols:
        ancestry[p] = ancestry[p].apply(lambda x: GetFloatFromPC(x), 1)
    return ancestry

def WriteGWAS(gwas, outpath, covars):
    # Ouptut header with command used
    f = open(outpath, "w")
    f.write("#" + " ".join(sys.argv) + "\n")
    f.write("# covars: %s\n"%",".join(covars))
    f.close()
    # Append gwas results
    gwas.to_csv(outpath, sep="\t", mode="a", index=False)

def Inverse_Quantile_Normalization(M):
    M = M.transpose()
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    Q = Q.transpose() 
    return Q

def NormalizeData(data, norm):
    # Add normalization quantile
    if norm == "quantile":
        data["phenotype"] = Inverse_Quantile_Normalization(data[["phenotype"]])
        return data

    # Add z-score normalization
    elif norm == "zscore":
        data["phenotype"]  = stats.zscore(data[["phenotype"]])
        return data

    else:
        ERROR("No normalization method specified")

def get_subregions(region, window_size):
    #region eg chr11:540543-678907
    #window_size default = 1MB
    chrom, rest = region.split(':')
    start, end = rest.split('-')
    start = int(start)
    end = int(end)
    num_windows = math.ceil((end - start) / window_size)
    windows = np.repeat(window_size, num_windows+1)
    windows[0] = start
    windows = np.cumsum(windows)
    if windows[-1] != end:
        windows[-1] = end
    subregions = [f'{chrom}:{p1}-{p2}' for p1,p2 in zip(windows[:-1], windows[1:])]
    return subregions

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--isbinary", help="True if phenotype is case control, False otherwise", type=bool, default=False)
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions", default=ANCESTRY_PRED_PATH)
    parser.add_argument("--region", help="chr:start-end to restrict to. Default is genome-wide", type=str)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--plot", help="Make a Manhattan plot", action="store_true")
    parser.add_argument("--norm", help="Normalize phenotype either quantile or zscore", type=str)
    parser.add_argument("--norm-by-sex",
                        help="Apply the normalization for each sex separately. Default: False",
                        action="store_true")
    parser.add_argument("--sample-call-rate", help="Apply minimum sample call rate QC", type=float, default=0.90)
    parser.add_argument("--variant-call-rate", help="Apply minimum variant call rate QC", type=float, default=0.90)
    parser.add_argument("--MAF", help="Apply minor allele frequency QC", type=float, default=0.01)
    parser.add_argument("--HWE", help="Apply HWE p-value cutoff QC", type=float, default=1e-100)
    parser.add_argument("--GQ", help="Apply minimun genotype score QC", type=int, default=20)
    args = parser.parse_args()

    # Set up paths
    if args.phenotype.endswith(".csv"):
        ptcovar_path = args.phenotype
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)

    # Check options
    if not CheckRegion(args.region):
        ERROR("Invalid region %s"%args.region)
    if args.num_pcs > 10:
        ERROR("Specify a maximum of 10 PCs")
    if args.norm_by_sex and args.norm is None:
        ERROR("Must specify --norm if using --norm-by-sex")

    # Get covarlist
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""]
    pt_covars = [item for item in args.ptcovars.split(",") if item != ""]
    covars = pcols + pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    data = pd.read_csv(ptcovar_path)
    #ancestry = LoadAncestry(args.ancestry_pred_path)
    #data = pd.merge(data, ancestry[["person_id"]+pcols], on=["person_id"])
    data["person_id"] = data["person_id"].apply(str)

    # Add normalization. If indicated, normalize for each sex separately.
    if args.norm_by_sex:
        # Separate the data into two smaller dataframes based on sex at birth.
        female_data = data[data['sex_at_birth_Male'] == 0].copy()
        male_data = data[data['sex_at_birth_Male'] == 1].copy()

        # Apply normalization on female and male dataframes separately.
        female_data =NormalizeData(data=female_data, norm=args.norm)
        male_data = NormalizeData(data=male_data, norm=args.norm)

        # Concatenate the female and male dataframes back into one
        # and sort the dataframe by original order.
        data = pd.concat([female_data, male_data])
        data = data.sort_index()
    else:
        # Apply normalization on the entire data.
        if args.norm is not None:
            data = NormalizeData(data=data, norm=args.norm)
        
    # Add shared covars
    sampfile = args.samples
    if sampfile.startswith("gs://"):
        sampfile = sampfile.split("/")[-1]
        if not os.path.isfile(sampfile):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(args.samples))
    samples = pd.read_csv(sampfile)
    samples["person_id"] = samples["person_id"].apply(str)
    
    data = pd.merge(data, samples)

    # Check we have all covars
    req_cols = ["phenotype"] + covars
    for item in req_cols:
        if item not in data.columns:
            ERROR("Required column %s not found"%item)


    # Set up GWAS method
    from hail_runner import HailRunner
    hl.init(default_reference = "GRCh38")
    window_size = 1000000
    subregions = get_subregions(args.region, window_size)
    print(subregions)
    all_subgwas = []
    for subregion in subregions:
        runner = HailRunner(data, isbinary=args.isbinary, region=subregion, covars=covars, sample_call_rate=args.sample_call_rate, variant_call_rate=args.variant_call_rate, MAF=args.MAF, HWE=args.HWE, GQ=args.GQ)
        # Run GWAS
        outpath = GetOutPath(args.phenotype, args.region, sampfile)
        runner.RunGWAS()
        if runner.gwas is not None:
            all_subgwas.append(runner.gwas)
    
    full_gwas = pd.concat(all_subgwas, axis=0).reset_index(drop=True)
    WriteGWAS(full_gwas, outpath+".tab", covars)

    # Plot Manhattan
    if args.plot:
        PlotManhattan(runner.gwas, outpath+".manhattan.png")
        PlotQQ(runner.gwas, outpath+".qq.png")

if __name__ == "__main__":
    main()
