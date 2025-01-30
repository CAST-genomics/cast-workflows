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
from cyvcf2 import VCF
import gzip
from io import StringIO
from collections import Counter

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

def GetCohort(samplefile):
    cohort = "ALL" #if samplefile == DEFAULTSAMPLECSV
    if samplefile != DEFAULTSAMPLECSV:
        cohort = os.path.basename(samplefile[:-4]) #chop off .csv at the end
    return cohort

def GetOutPath(phenotype, method, region, samplefile, outdir):
    cohort = GetCohort(samplefile)
    outprefix = "%s_%s_%s"%(phenotype, method, cohort)
    if region is not None:
        outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    return os.path.join(outdir, outprefix + ".gwas")

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
        if line.startswith("chrom"):
            # Header line
            continue
        chrom, start, gene = line.strip().split(",")[:3]
        annotations.append((str(chrom), str(start), str(gene)))
    return annotations

def get_alleles(record):
    ref_allele = record.REF
    alt_alleles = record.ALT
    info_fields = record.INFO
    ru = info_fields["RU"]
    if ru is None:
        print("Error: Repeat Unit (RU) not provided in the info field")
        return
    ru_len = len(ru)
    #print("ru_len, ru, len ref allele", ru_len, ru, len(ref_allele))
    len_ref_allele = len(ref_allele)
    len_alt_alleles = [len(alt_allele) for alt_allele in alt_alleles]
    return len_ref_allele, len_alt_alleles, ru_len

def get_rc_from_allele(allele, ref_allele, alt_alleles, ru_len):
    if allele == 0:
        ref_rc = int(ref_allele/ru_len)
        #print("ref_rc ", ref_rc)
        return ref_rc
    else:
        rc = int(alt_alleles[allele-1]/ru_len)
        #print("call {} allele {}".format(allele, rc))
        return rc

def set_genotypes(data, annotations, cohort, samples, phenotype, tr_vcf, outdir):
    print("Reading genotypes from the vcf file")
    # Extract the chrom from the tr_vcf file
    words = re.split("_|\.", tr_vcf)
    vcf_chrom = [word for word in words if word.startswith("chr")][0]
    # Plot phenotype histogram
    plot_histogram(data["phenotype"], phenotype,
                    os.path.join(outdir,
                        "{}_histogram_after_norm.png".format(phenotype)))
    samples_str = [str(sample) for sample in samples]
    vcf = VCF(tr_vcf, samples = samples_str)
    vcf_samples = vcf.samples

    # Focus on a few loci of interest
    for chrom, start, gene in annotations:
        all_alleles = []
        empty_calls, no_calls = 0, 0
        if chrom != vcf_chrom:
            continue
        variant = list(vcf("{}:{}-{}".format(chrom, int(start)-10, int(start) + 10)))[0]
        print("Processing annotation {} {}:{}".format(gene, chrom, start))
        data.loc[:, [gene]] = np.nan
        samples_with_calls = set()
        print("Reading ref and alt alleles")
        ref_allele_len, alt_alleles_len, ru_len = get_alleles(variant)
        print("Reading calls")
        counter = 0
        for i, sample_call in enumerate(variant.genotypes):
            counter += 1
            sample = vcf_samples[i]
            if counter % 10000 == 0:
                print("reading call ", counter)
                break
            rc_1 = get_rc_from_allele(int(sample_call[0]), ref_allele_len, alt_alleles_len, ru_len)
            rc_2 = get_rc_from_allele(int(sample_call[1]), ref_allele_len, alt_alleles_len, ru_len)
            rc = (rc_1 + rc_2) / 2.0
            all_alleles.extend([rc_1, rc_2])

            samples_with_calls.add(sample)
            data.loc[data["person_id"]==sample, gene] = rc
        print("Alleles count for the 4 most common alleles: ", Counter(all_alleles).most_common(4))
        # Plot alleles
        if len(set(all_alleles)) > 1:
            print("Plotting alleles histogram")
            print("gene: ", gene)
            print("phenotype: ", phenotype)
            # Polymorphic vntr in the imputed set. Otherwise, if it's non-polymorphic, it'll get an error.
            plot_histogram(all_alleles, gene, os.path.join(outdir, "{}_alleles_{}_{}.png".format(gene, phenotype, cohort)))
        if no_calls + empty_calls > 0:
            print("Skipping {} empty calls and {} no calls for {} on vcf".format(
                    empty_calls, no_calls, gene))
    return data

def print_stats(data, locus):
    print("For phenotype, mean {:.4f} and sd {:.4f}".format(
        data["phenotype"].mean(),
        data["phenotype"].std()))
    print(data["phenotype"].describe())
    print("For genotype, mean {:.4f} and sd {:.4f}".format(
        data[locus].mean(),
        data[locus].std()))
    print(data[locus].describe())
    min_height = data["phenotype"].min()
    max_height = data["phenotype"].max()
    samples_w_min_height = data[data["phenotype"]==min_height]
    print("samples_w_min_height: ", len(samples_w_min_height))
    samples_w_max_height = data[data["phenotype"]==max_height]
    print("samples_w_max_height: ", len(samples_w_max_height))
    
    print("median genotype at samples_w_min_height: ", samples_w_min_height[locus].median())
    print("median genotype at samples_w_max_height: ", samples_w_max_height[locus].median())

def load_snp_gwas_df(filename):
    snp_data = pd.read_csv(filename)
    print("Num snp variants loaded: ", len(snp_data))
    snp_data["POS"] = snp_data["POS"].astype(int)
    snp_data = snp_data.rename(columns={"POS": "pos",
                     "Pvalue_log10": "-log10pvalue",
                     "CHR": "chrom"})
    return snp_data

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--method", help="GWAS method. Options: %s"",".join(GWAS_METHODS), type=str, default="hail")
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions", default=ANCESTRY_PRED_PATH)
    parser.add_argument("--region", help="chr:start-end to restrict to. Default is genome-wide", type=str)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age",
                            type=str, default="age")
    parser.add_argument("--sharedcovars",
                help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male",
                            type=str, default="sex_at_birth_Male")
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes. Required if running associaTR", type=str)
    parser.add_argument("--is-imputed", help="Indicate if the tr-vcf file is imputed", action="store_true")
    parser.add_argument("--annotations", help="CSV file with annotations for plotting" + \
                                              "each line should include the " + \
                                              "chromosome, start coordinate and label(gene).",
                                              type=str)
    parser.add_argument("--snp-gwas-file", help="File where SNP gwas dataframe was stored from the "+ \
                                                "all by all dataset")
    parser.add_argument("--plot", help="Make a Manhattan plot", action="store_true")
    parser.add_argument("--norm", help="Normalize phenotype either quantile or zscore", type=str)
    parser.add_argument("--norm-by-sex",
                        help="Apply the normalization for each sex separately. Default: False",
                        action="store_true")
    parser.add_argument("--outdir", help="Path to the output directory", type=str, required=True)
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

    # Create output dir if does not exist
    if os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

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

    # Load SNP gwas
    snp_gwas = None
    if args.snp_gwas_file is not None:
        snp_gwas = load_snp_gwas_df(args.snp_gwas_file)

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

    # Add shared covars
    sampfile = args.samples
    if sampfile.startswith("gs://"):
        sampfile = sampfile.split("/")[-1]
        if not os.path.isfile(sampfile):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(args.samples))
    samples = pd.read_csv(sampfile)
    samples["person_id"] = samples["person_id"].apply(str)

    # Set cohort name which is later used to name output files
    cohort = GetCohort(sampfile)
    
    # Plot genotype-phenotype plot and allele histogram
    if args.method == "associaTR":
        if args.annotations is not None:
            # Get annotations of specific TRs for plotting
            annotations = read_annotations(args.annotations)
            # This is necessary for genotype-phenotype plot
            # To run a quick check, randomly subsample the samples
            data = set_genotypes(data=data,
                                 annotations=annotations,
                                 cohort=cohort,
                                 samples=samples["person_id"],
                                 phenotype=args.phenotype,
                                 tr_vcf=args.tr_vcf,
                                 outdir=args.outdir)
    
    data = pd.merge(data, samples)

    # Check we have all covars
    print("Check all covars are present")
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
    outpath = GetOutPath(args.phenotype, args.method, args.region, sampfile, args.outdir)
    runner.RunGWAS()
    WriteGWAS(runner.gwas, outpath+".tab", covars)

    # Plot Manhattan
    if args.plot:
        if args.method == "associaTR":
            annotate = True
            p_value_threshold = -np.log10(5*10**-8)
            if args.annotations is not None:
                #print("plotting genotype phenotype for annotations ", annotations)
                for chrom, pos, gene in annotations:
                    #print("chrom, pos and gene", chrom, pos, gene)
                    plot_genotype_phenotype(data=data,
                        genotype=gene,
                        gwas=runner.gwas,
                        chrom=chrom,
                        pos=pos,
                        phenotype="phenotype",
                        phenotype_label=args.phenotype,
                        outpath=os.path.join(args.outdir,
                                "{}_genotype_{}_{}.png".format(
                                    gene, args.phenotype, cohort)))
        else:
            # no text annotation on manhattan plot for hail runner
            annotate = False
            p_value_threshold = -np.log10(5*10**-8)

        gwas = runner.gwas
        # Combine two dataframes into one, with source indicated in a new column.
        columns = ["chrom", "pos", "-log10pvalue"]
        if snp_gwas is not None:
            columns = ["chrom", "pos", "-log10pvalue"]
            gwas = gwas[columns]
            gwas["variant"] = "VNTR"
            snp_gwas = snp_gwas[columns]
            snp_gwas["variant"] = "SNP"
            print("gwas shape", gwas.shape)
            print("snp_gwas shape", snp_gwas.shape)
            gwas = pd.concat([
                        snp_gwas[snp_gwas.index.isin(gwas.index) == False],
                        gwas])
            print("combined gwas shape", gwas.shape)
            columns = ["variant"] + columns
            hue = "variant"
        else:
            hue = "chrom"


        # Store significant hits
        significant_hits = gwas[gwas["-log10pvalue"] > p_value_threshold]
        significant_hits.sort_values(by="pos", ignore_index=True, inplace=True)
        significant_hits = significant_hits[columns]
        significant_hits.to_csv(outpath + "_significant_hits.csv", index=False)

        PlotManhattan(gwas, outpath+".manhattan.png",
                      hue=hue,
                      annotations=annotations,
                      annotate=annotate,
                      p_value_threshold=p_value_threshold)
        PlotQQ(runner.gwas, outpath+".qq.png")

if __name__ == "__main__":
    main()
