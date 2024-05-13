#!/usr/bin/env python3

"""
Run PheWAS on All of Us
Example: 
./aou_phewas.py --tr-vcf CBL_test.filtered.sorted.vcf.gz --region chr11:119206290-119206323 --out CBL_phewas
"""

import argparse
import pandas as pd
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from statsmodels.regression.linear_model import OLS
import statsmodels.formula.api as smf
import statsmodels.api as sm
import sys
import os
import numpy as np
import scipy.stats as stats
from utils import MSG, ERROR

ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")
MANIFESTFILE ="https://github.com/CAST-genomics/cast-workflows/raw/nichole/concept_id_delete_Kelvin/gwas/aou/phenotypes_manifest.csv"
#MANIFESTFILE = "https://github.com/CAST-genomics/cast-workflows/raw/main/gwas/aou/phenotypes_manifest.csv"

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

def Inverse_Quantile_Normalization(M):
    M = M.transpose()
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    Q = Q.transpose() 
    return Q

# def NormalizeData_Quantile(data):
#     # Add normalization quantile
#     data["phenotype"] = Inverse_Quantile_Normalization(data[["phenotype"]])
#     data["age"] = Inverse_Quantile_Normalization(data[["age"]])
#     return data


# def NormalizeData_Zscore(data):
#     # Add z-score normalization
#     data["phenotype"]  = stats.zscore(data[["phenotype"]])
#     data["age"]  = stats.zscore(data[["phenotype"]])
#     return data

 
def main():
    parser = argparse.ArgumentParser(__doc__)
    # Specifying the phenotypes
    parser.add_argument("--manifest", help="Manifest file with phenotype info", type=str, default=MANIFESTFILE)
    # Specifying the genotypes
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes.", type=str)
    parser.add_argument("--region", help="chr:start-end of target variant", type=str, required=True)
    # Options for shared covariates
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions", default=ANCESTRY_PRED_PATH)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--quantile", help="Quantile normalize phenotype", action="store_true", default=False)   
    parser.add_argument("--zscore", help="Zscore normalize phenotype", action="store_true", default=False)
    # Output options
    parser.add_argument("--out", help="Name of output file", type=str, default="stdout")
    args = parser.parse_args()

    # Set up output file
    if args.out == "stdout":
        outf = sys.stdout
    else:
    	outf = open(args.out, "w")

    # Set up shared covariates
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""] + pcols
    data = LoadAncestry(args.ancestry_pred_path)
    data["person_id"] = data["person_id"].apply(str)


    # Load TR genotypes for the target locus to a df and merge with data
    # Genotype should be in a column labeled "genotype"
    # TODO , check in with Melissa, what genotype info to use, so far use allele_sum

    # Load TR genotypes
    invcf = utils.LoadSingleReader(args.tr_vcf, checkgz=True)
    samples = invcf.samples
    region = invcf(args.region)
    nrecords = 0
    for record in region:
        trrecord = trh.HarmonizeRecord(trh.VcfTypes["hipstr"], record)
        afreqs = trrecord.GetAlleleFreqs()
        genotypes = trrecord.GetLengthGenotypes()
        allele_sum = [sum(item)for item in genotypes]
        trdf = pd.DataFrame({"person_id": samples, "genotype": allele_sum})
        nrecords += 1
    if nrecords == 0:
        ERROR("No matching TR records found")
    if nrecords > 1:
        ERROR("Multiple matching TR records found")

    #merge genotype with data
    data = pd.merge(data, trdf, on=["person_id"])
    
    # Process the phenotypes from manifest file one at a time
    manifest = pd.read_csv(args.manifest)
    for index, row in manifest.iterrows():
        phenotype = row["phenotype"]
    	# Load phenotype data
        ptfile = os.path.join(os.environ["WORKSPACE_BUCKET"],"phenotypes", row["phenotype_file"])
        ptdata = pd.read_csv(ptfile)
        ptdata["person_id"] = ptdata["person_id"].apply(str)
        ptcovars = [item for item in ptdata.columns if item != phenotype]
    	# Merge with genotypes. add intercept
        ptdata = pd.merge(data, ptdata, on=["person_id"])
        #add normalization of phenotype and covariates
        if args.quantile:
            ptdata["phenotype"] = Inverse_Quantile_Normalization(ptdata[["phenotype"]])
            ptdata["age"] = Inverse_Quantile_Normalization(ptdata[["age"]])
            for i in range(1,16):
                pc_name = f'PC_{i}'
                ptdata[pc_name] = Inverse_Quantile_Normalization(ptdata[[pc_name]])
            
            MSG("Quantile normalizing phenotype,age and pcs")

        if args.zscore:
            ptdata["phenotype"]  = stats.zscore(ptdata[["phenotype"]])
            ptdata["age"]  = stats.zscore(ptdata[["phenotype"]])
            for i in range(1,16):
                pc_name = f'PC_{i}'
                ptdata[pc_name] = stats.zscore(ptdata[[pc_name]])

            MSG("Zscore normalizing phenotype,age and pcs")
        ptdata["intercept"] = 1
        # Regressionls
        covars = ["intercept"] + shared_covars + ptcovars 
    	# TODO, redo - need to put a flag in manifest to know if something
        #if logistic:
            #model = smf.logit(ptdata["phenotype"], ptdata[["genotype"]+covars].astype(float))

        model = sm.OLS(ptdata["phenotype"], ptdata[["genotype"]+covars].astype(float))
            
        reg_result = model.fit()
        pval = reg_result.pvalues[0]
        coef = reg_result.params[0]
        se = reg_result.bse[0]
        outf.write("\t".join([phenotype, str(pval), str(coef), str(se)])+"\n")
    outf.close()

if __name__ == "__main__":
    main()
