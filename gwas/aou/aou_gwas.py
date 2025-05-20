#!/usr/bin/env python3

"""
Run GWAS on All of Us

Example:
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
"""

import argparse
import os
import pandas as pd
import re
import sys
from utils import MSG, ERROR
import numpy as np
from gwas_plotter import PlotManhattan, PlotQQ, plot_histogram
from plot_genotype_phenotype import plot_genotype_phenotype
from gwas_utils import SAMPLEFILE, GetPTCovarPath, GetCohort, CheckRegion, NormalizeData, \
                        read_annotations, set_genotypes, load_snp_gwas_df, load_gwas_tab

GWAS_METHODS = ["hail", "associaTR"]
ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"


def GetOutPath(phenotype, method, region, samplefile, outdir):
    cohort = GetCohort(samplefile)
    outprefix = "%s_%s_%s"%(phenotype, method, cohort)
    if region is not None:
        outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    return os.path.join(outdir, outprefix + ".gwas")

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
        ancestry[p] = ancestry[p].str.replace("[","").str.replace("]","").astype(float)
    return ancestry

def WriteGWAS(gwas, outpath, covars):
    # Ouptut header with command used
    f = open(outpath, "w")
    f.write("#" + " ".join(sys.argv) + "\n")
    f.write("# covars: %s\n"%",".join(covars))
    f.close()
    # Append gwas results
    gwas.to_csv(outpath, sep="\t", mode="a", index=False)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--method", help="GWAS method. Options: %s"",".join(GWAS_METHODS), type=str, default="hail")
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions", default=ANCESTRY_PRED_PATH)
    parser.add_argument("--region", help="chr:start-end to restrict to. Default is genome-wide", type=str)
    parser.add_argument("--binary", help="Plot the genotype-phenotype for binary phenotypes", action="store_true")
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age",
                            type=str, default="age")
    parser.add_argument("--sharedcovars",
                help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male",
                            type=str, default="sex_at_birth_Male")
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes. Required if running associaTR", type=str)
    parser.add_argument("--annotations", help="CSV file with annotations for plotting" + \
                                              "each line should include the " + \
                                              "chromosome, start coordinate and label(gene).",
                                              type=str)
    parser.add_argument("--snp-gwas-file", help="File where SNP gwas dataframe was stored from the "+ \
                                                "all by all dataset")
    parser.add_argument("--vntr-gwas-file", help="File where VNTR gwas dataframe was stored from the "+ \
                                                "plink run, usually needed for binary traits.")
    parser.add_argument("--vntr-as-covariate", help="VNTR (chrom_start_gene) to be considered as "+ \
                                                "a covariate for a conditional snp gwas.")
    parser.add_argument("--plot", help="Make a Manhattan plot", action="store_true")
    parser.add_argument("--violin", help="Make the genotype_phenotype plot a violinplot instead of boxplot.",
                            action="store_true")
    parser.add_argument("--merge-bins", help="Merge bins in the genotype_phenotype plot.",
                            action="store_true")
    parser.add_argument("--plot-genotype-phenotype", help="Plot the genotype phenotype joint plot. " + \
                        "This may take about half an hour for each single annotation point in the chromosome.",
                        action="store_true")
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
    ptcovar_path = "phenocovars/{}_phenocovar.csv".format(args.phenotype)
    if os.path.exists(ptcovar_path):
        pass
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)

    if args.violin and args.plot_genotype_phenotype is False:
        ERROR("--violin must be used when --plot-genotype-phenotype is set.")
    if args.merge_bins and args.plot_genotype_phenotype is False:
        ERROR("--merge-bins must be used when --plot-genotype-phenotype is set.")
    if args.plot_genotype_phenotype and args.annotations is None:
        ERROR("--annotations must be provided when --plot-genotype-phenotype is set.")
        
    # Create output dir if does not exist
    if not os.path.exists(args.outdir):
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
    if args.vntr_as_covariate and args.method != "hail":
        ERROR("--vntr-as-covariate could only be used with SNP-gwas. Use hail for --method.")

    # Load SNP gwas from all by all query or VNTR gwas from a previous run
    snp_gwas = None
    if args.snp_gwas_file:
        snp_gwas = load_snp_gwas_df(args.snp_gwas_file)
    if args.vntr_gwas_file:
        vntr_gwas = load_gwas_tab(args.vntr_gwas_file)

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
    
    # Get annotations of specific TRs for plotting
    annotations = read_annotations(args.annotations)

    # Extract the chrom from the tr_vcf file for the manhattan plot
    words = re.split("_|\.", args.tr_vcf)
    vcf_chrom = [word for word in words if word.startswith("chr")][0]

    # Create a more natural phenotype label
    phenotype_label = args.phenotype.lower().replace("_", " ")
    if phenotype_label == "red blood cell distribution width":
        phenotype_label = "RDW"

    # Get repeat counts from the tr-vcf file and plot alleles histogram
    if args.method == "associaTR" \
            and args.plot_genotype_phenotype:
        # Setting genotypes is necessary for genotype-phenotype plot.
        # The first time computing it for each phenotype-VNTR locus pair
        # takes some time (about an hour). Consecutive times are much faster
        # as genotype data is automatically stored the first time computed and loaded later.
        data, annotations = set_genotypes(data=data,
                         annotations=annotations,
                         cohort=cohort,
                         samples=samples["person_id"],
                         tr_vcf=args.tr_vcf,
                         outdir=args.outdir)
        # Plot phenotype histogram
        suffix = "no_norm"
        if args.norm:
            suffix = "after_norm"
        plot_histogram(data["phenotype"], phenotype_label,
                    os.path.join(args.outdir,
                        "{}_histogram_{}.png".format(args.phenotype, suffix)))
    # Using vntr as covariate for a conditional snp gwas.
    if args.vntr_as_covariate:
        annotations = [args.vntr_as_covariate.split("_")]
        gene = annotations[0][2]
        g_data, annotations = set_genotypes(data=data,
                         annotations=annotations,
                         cohort=cohort,
                         samples=samples["person_id"],
                         tr_vcf=args.tr_vcf,
                         outdir=args.outdir)
        covars = covars + [gene]
        data[gene] = g_data[gene]

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
    if args.vntr_gwas_file:
        gwas = vntr_gwas.loc[vntr_gwas["chrom"]==vcf_chrom]
    else:   
        runner.RunGWAS()
        WriteGWAS(runner.gwas, outpath+".tab", covars)
        gwas = runner.gwas

    # Plot Manhattan
    if args.plot:
        if args.method == "associaTR":
            annotate = True
            p_value_threshold = -np.log10(5*10**-8)
            if args.plot_genotype_phenotype:
                for chrom, pos, gene in annotations:
                    plot_genotype_phenotype(data=data,
                        genotype=gene,
                        phenotype="phenotype",
                        phenotype_label=phenotype_label,
                        binary=args.binary,
                        violin=args.violin,
                        merge_bins=args.merge_bins,
                        outpath=os.path.join(args.outdir,
                                "{}_genotype_{}_{}.png".format(
                                    gene, args.phenotype, cohort)))
                    if args.binary:
                        effect_column = "OR"
                    else:
                        effect_column = "beta"
                    if len(gwas) == 1:
                        effect_size = gwas.iloc[0][effect_column]
                        p_val = gwas.iloc[0]["-log10pvalue"]
                    else: # To avoid warning on single element series
                        effect_size = gwas.loc[(gwas["chrom"] == chrom) &\
                                        (gwas["pos"] > int(pos) - 1) &\
                                        (gwas["pos"] < int(pos) + 1), effect_column].item()
                        p_val = gwas.loc[(gwas["chrom"] == chrom) &\
                                        (gwas["pos"] > int(pos) - 1) &\
                                        (gwas["pos"] < int(pos) + 1)]["-log10pvalue"].item()
                    print("Effect ({}) for {} is {} and -log10p is {}".format(
                            effect_column, gene, effect_size, p_val))
        else:
            # no text annotation on manhattan plot for hail runner
            annotate = False
            p_value_threshold = -np.log10(5*10**-8)

        if args.region:
            chrom, start_end = args.region.split(":")
            start, end = start_end.split("-")
            start, end = int(start), int(end)
            gwas = gwas[(gwas["chrom"] == chrom) & \
                        (gwas["pos"] >= start) & \
                        (gwas["pos"] <= end)]
        # Combine two dataframes into one, with source indicated in a new column.
        columns = ["chrom", "pos", "-log10pvalue"]
        if snp_gwas is not None:
            columns = ["chrom", "pos", "-log10pvalue"]
            gwas = gwas[columns]
            gwas.loc[:,"variant"] = "VNTR"
            snp_gwas = snp_gwas[columns]
            snp_gwas.loc[:,"variant"] = "SNP"
            print("gwas shape", gwas.shape)
            print("snp_gwas shape", snp_gwas.shape)
            gwas = pd.concat([snp_gwas, gwas], ignore_index=True)
            print("combined gwas shape", gwas.shape)
            columns = ["variant"] + columns
            hue = "variant"
        else:
            hue = "chrom"


        # Store significant hits
        significant_hits = gwas.loc[gwas["-log10pvalue"] > p_value_threshold]
        significant_hits = significant_hits.sort_values(by="pos", ignore_index=True)
        significant_hits = significant_hits[columns]
        significant_hits.to_csv(outpath + "_significant_hits.csv", index=False)

        PlotManhattan(gwas, outpath+".manhattan.png",
                      hue=hue,
                      phenotype=phenotype_label,
                      annotations=annotations,
                      annotate=annotate,
                      p_value_threshold=p_value_threshold)
        PlotQQ(gwas, outpath+".qq.png")

if __name__ == "__main__":
    main()
