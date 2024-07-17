#!/usr/bin/env python3

"""
Run GWAS on All of Us

Example:
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
"""

import argparse
from gwas_plotter import PlotManhattan, PlotQQ, plot_histogram, plot_genotype_phenotype
import os
import pandas as pd
import re
import sys
from utils import MSG, ERROR
import numpy as np
import scipy.stats as stats
import warnings
import re
import vcf
import gzip
from io import StringIO

GWAS_METHODS = ["hail", "associaTR"]
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

def GetOutPath(phenotype, method, region, samplefile):
    cohort = "ALL" #if samplefile == DEFAULTSAMPLECSV
    if samplefile != DEFAULTSAMPLECSV:
        cohort = samplefile[:-4] #chop off .csv at the end
    outprefix = "%s_%s_%s"%(phenotype, method, cohort)
    if region is not None:
        outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
    return "outputs/" + outprefix + ".gwas"

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

def read_annotations(filename):
    annotations = []
    with open(filename, "r") as annotations_file:
        lines = annotations_file.readlines()
    for line in lines:
        chrom, start, gene = line.strip().split(",")
        annotations.append((str(chrom), str(start), str(gene)))
    return annotations

def get_rc_from_allele_record(record, allele):
    ref_allele = record.REF
    alt_alleles = record.ALT.split(",")
    ru = record.INFO.split(";")[2].replace("RU=", "")
    ru_len = len(ru)
    #print("ru_len, ru, len ref allele", ru_len, ru, len(ref_allele))
    if allele == 0:
        ref_rc = int(len(ref_allele)/ru_len)
        #print("ref_rc ", ref_rc)
        return ref_rc
    else:
        rc = int(len(alt_alleles[allele-1])/ru_len)
        #print("call {} allele {}".format(allele, rc))
        return rc

def get_mean_rc_from_call_record(record, call, sep="/"):
    allele_1, allele_2 = call.split(sep)
    rc_1 = get_rc_from_allele_record(record, int(allele_1))
    rc_2 = get_rc_from_allele_record(record, int(allele_2))
    return (rc_1 + rc_2) / 2

def get_rc_from_allele(vcf_df_row, allele):
    ref_allele = vcf_df_row["REF"].array[0]
    alt_alleles = vcf_df_row["ALT"].array[0].split(",")
    ru = vcf_df_row["INFO"].array[0].split(";")[2].replace("RU=", "")
    ru_len = len(ru)
    #print("ru_len, ru, len ref allele", ru_len, ru, len(ref_allele))
    if allele == 0:
        ref_rc = int(len(ref_allele)/ru_len)
        #print("ref_rc ", ref_rc)
        return ref_rc
    else:
        rc = int(len(alt_alleles[allele-1])/ru_len)
        #print("call {} allele {}".format(allele, rc))
        return rc

def get_individual_rcs_from_call(vcf_df_row, call, sep="/"):
    allele_1, allele_2 = call.split(sep)
    rc_1 = get_rc_from_allele(vcf_df_row, int(allele_1))
    rc_2 = get_rc_from_allele(vcf_df_row, int(allele_2))
    return rc_1, rc_2

def get_overall_rc_from_call(vcf_df_row, call, sep="/"):
    allele_1, allele_2 = call.split(sep)
    rc_1 = get_rc_from_allele(vcf_df_row, int(allele_1))
    rc_2 = get_rc_from_allele(vcf_df_row, int(allele_2))
    return (rc_1 + rc_2)/2.0


def set_genotypes(data, args, annotations):
    #print("is imputed: ", args.is_imputed)
    imputed = args.is_imputed
    # Plot phenotype histogram
    plot_histogram(data["phenotype"], "outputs/{}_histogram_after_norm.png".format(args.phenotype))
    # Read input VCF file into a dataframe
    lines = None

    if args.tr_vcf.endswith("gz"):
        # It is a compressed vcf file
        with gzip.open(args.tr_vcf, "rt") as vcf_file:
            lines = "\n".join(vcf_file.readlines())
    elif args.tr_vcf.endswith("vcf"):
        # It is an uncompressed vcf file
        with open(args.tr_vcf, "rb") as vcf_file:
            lines = "\n".join(vcf_file.readlines())
    else:
        print("Error: Cannot recognize tr-vcf file format. Should be either a vcf file or a vcf.gz file")
        exit(1)
    vcf_df = pd.read_csv(StringIO(re.sub("#CHROM", "CHROM", lines)), sep="\t", comment='#')

    # Set up an output vcf file for normalized values
    #normalized_file = open("normalized.vcf", "w+")
    #for line in lines:
    #    if line.startswith("#"):
    #        normalized_file.write(line + "\n")
    # Focus on a few loci of interest
    for chrom, start, gene in annotations:
        all_alleles = []
        empty_calls, no_calls = 0, 0

        locus_calls = vcf_df[(vcf_df["CHROM"] == chrom) & \
                             (vcf_df["POS"] == int(start))
                             ]
        data.loc[:, [gene]] = np.nan
        samples_with_calls = set()
        shared_columns = []
        for column in vcf_df.columns:
            if column.isnumeric():
                # Corresponds to a sample id
                if not imputed:
                    #print("Column: ", column)
                    #print("Locus call", locus_calls[column].to_string())
                    if len(locus_calls[column].array) == 0:
                        # An error causing an empty value in the vcf file
                        empty_calls += 1
                        continue
                    call = locus_calls[column].array[0]
                    call = call.split(":")[0]
                    if call == ".":
                        # Equal to no call
                        no_calls += 1
                        continue
                    rc = get_overall_rc_from_call(locus_calls, call, sep="/")
                    # Get individual alleles for plotting
                    rc_1, rc_2 = get_individual_rcs_from_call(locus_calls, call, sep="/")
                    all_alleles.extend([rc_1, rc_2])
                else: # imputed
                    if len(locus_calls[column]) == 0:
                        # Equal to no call
                        continue
                    call = locus_calls[column].array[0]
                    call = call.split(":")[0]
                    #print("Locus call after split ", call)
                    rc = get_overall_rc_from_call(locus_calls, call, sep="|")
                    # Get individual alleles for plotting
                    rc_1, rc_2 = get_individual_rcs_from_call(locus_calls, call, sep="|")
                    all_alleles.extend([rc_1, rc_2])


                # Remove outliers for CACNA1C
                #if gene == "CACNA1C" and (rc < 2*180 or rc > 2*205):
                #    continue
                samples_with_calls.add(column)
                data.loc[data["person_id"]==column, gene] = rc
            else:
                shared_columns.append(column)
        #data = data.dropna(subset=[gene, "phenotype"])
        
        # Plot individual alleles for ACAN
        if gene == "ACAN" and args.phenotype == "height":
            plot_histogram(all_alleles, "outputs/ACAN_alleles_height.png")
        if no_calls + empty_calls > 0:
            print("Skipping {} empty calls and {} no calls for {} on vcf".format(
                    empty_calls, no_calls, gene))
        # Normalize genotypes
        #data.loc[:, [gene]] = stats.zscore(data[gene])
        #data.loc[:, [gene]] = Inverse_Quantile_Normalization(data[[gene]])
        
        # Write output normalized file
        #locus_normalized = data.loc[:, [gene]]
        #normalized_file.write(shared_columns + )


        #print("len of data after removing samples with no calls or no phenotypes {}: {}".format(
        #        args.phenotype, len(data)))
        #print(plotted_data.head())
        # Plotting genotype phenotype later when we have the effect sizes.
        #plot_genotype_phenotype(data=data,
        #            genotype=gene,
        #            phenotype="phenotype",
        #            outpath="outputs/{}_genotype_{}.png".format(gene, args.phenotype))
    return data


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--method", help="GWAS method. Options: %s"",".join(GWAS_METHODS), type=str, default="hail")
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions", default=ANCESTRY_PRED_PATH)
    parser.add_argument("--region", help="chr:start-end to restrict to. Default is genome-wide", type=str)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes. Required if running associaTR", type=str)
    parser.add_argument("--is-imputed", help="Indicate if the tr-vcf file is imputed", action="store_true")
    parser.add_argument("--annotations", help="CSV file with annotations for plotting" + \
                                              "each line should include the " + \
                                              "chromosome, start coordinate and label(gene).",
                                              type=str)
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
    ptcovar_path = "{}_phenocovar.csv".format(args.phenotype)
    #if args.phenotype.endswith(".csv"):
    if os.path.exists(ptcovar_path):
        pass
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)

    # Check options
    if args.method not in GWAS_METHODS:
        ERROR("--method must be one of: %s"%",".join(GWAS_METHODS))
    if not CheckRegion(args.region):
        ERROR("Invalid region %s"%args.region)
    if args.num_pcs > 10:
        ERROR("Specify a maximum of 10 PCs")
    if args.method == "associaTR" and args.tr_vcf is None:
        ERROR("Must specify --tr-vcf for associaTR")
    if args.norm_by_sex and args.norm is None:
        ERROR("Must specify --norm if using --norm-by-sex")

    # Get covarlist
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""]
    pt_covars = [item for item in args.ptcovars.split(",") if item != ""]
    covars = pcols + pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    data = pd.read_csv(ptcovar_path)
    ancestry = LoadAncestry(args.ancestry_pred_path)
    data = pd.merge(data, ancestry[["person_id"]+pcols], on=["person_id"])
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
    if args.method == "associaTR":
        # Get annotations of specific TRs for plotting
        annotations = read_annotations(args.annotations)
        data = set_genotypes(data, args, annotations)

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
    if args.method == "hail":
        from hail_runner import HailRunner
        runner = HailRunner(data, region=args.region, covars=covars, sample_call_rate=args.sample_call_rate, variant_call_rate=args.variant_call_rate, MAF=args.MAF, HWE=args.HWE, GQ=args.GQ)
    elif args.method == "associaTR":
        from associatr_runner import AssociaTRRunner
        runner = AssociaTRRunner(data, args.tr_vcf, region=args.region, covars=covars)
    else:
        ERROR("GWAS method %s not implemented" % args.method)

    # Run GWAS
    outpath = GetOutPath(args.phenotype, args.method, args.region, sampfile)
    runner.RunGWAS()
    WriteGWAS(runner.gwas, outpath+".tab", covars)
    #runner.gwas = pd.read_csv(outpath+".tab", sep="\t")

    # Plot Manhattan
    if args.plot:
        if args.method == "associaTR":
            annotate = True
            p_value_threshold = -np.log10(5*10**-3)
            for chrom, pos, gene in annotations:
                plot_genotype_phenotype(data=data,
                        genotype=gene,
                        gwas=runner.gwas,
                        chrom=chrom,
                        pos=pos,
                        phenotype="phenotype",
                        outpath="outputs/{}_genotype_{}.png".format(gene, args.phenotype))
        else:
            # no text annotation on manhattan plot for hail runner
            annotate = False
            p_value_threshold = -np.log10(5*10**-8)
            #print("GWAS index and column")
            #print(runner.gwas.index)
            #print(runner.gwas.columns)
            #added_points = {"chrom": ["chr15", "chr15"],
            #                "pos": [88855424,88857434],
            #                "p_value": [0.1, 0.01],
            #                }
            #runner.gwas = runner.gwas.append(added_points, ignore_index=True)
            #runner.gwas = runner.gwas.sort(['chrom', 'pos']).reset_index(drop=True)


        PlotManhattan(runner.gwas, outpath+".manhattan.png",
                      annotate=annotate,
                      p_value_threshold=p_value_threshold,
                      extra_points=[
                          ("chr15", 88855424, 1, "ACAN_v_s"),
                          ("chr15", 88857434, 1, "ACAN_v_e")
                          ])
        PlotQQ(runner.gwas, outpath+".qq.png")

if __name__ == "__main__":
    main()
