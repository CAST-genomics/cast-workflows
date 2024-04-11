#!/usr/bin/env python3

"""
Run PheWAS on All of Us

"""

import argparse
import pandas as pd
from statsmodels.regression.linear_model import OLS
import sys
import os

ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")
MANIFESTFILE = "https://github.com/CAST-genomics/cast-workflows/raw/main/gwas/aou/phenotypes_manifest.csv"

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
    # Output options
    parser.add_argument("--out", help="Name of output file", type=str, default="stdout")
    parser.add_argument("--logistic", help="Apply logistic regression",action="store_true", default=False)
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
    # TODO , try platelet count, which file to use??
    #genotype = pd.read_csv(args.tr_vcf)
    #print(genotype.head())

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
        ptdata["intercept"] = 1
    	# Regression
        #Nichole added
        #intercept = ptdata["intercept"]
        #covars = intercept + shared_covars + ptcovars
        print(ptdata.head())

    	# TODO - need to put a flag in manifest to know if something
        #if args.logistic:
            #add logistic regression
        #else:
          #add linear regresssion
    	# is case/control or quantitative. if case/control, use
    	# logistic regression instead
        model = OLS(ptdata["phenotype"], ptdata[["genotype"]+covars])
        reg_result = model.fit()
        pval = reg_result.pvalues[0]
        coef = reg_result.params[0]
        se = reg_result.bse[0]
        outf.write("\t".join([phenotype, str(pval), str(coef), str(se)])+"\n")
    outf.close()

if __name__ == "__main__":
    main()
