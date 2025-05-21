
"""
Plot genotype-phenotype plots and odds ratio plots depending on the phenotype being binary or continuous.
These plots are currently automatically created when running GWAS (aou_gwas.py), but can be generated separately as well.

Example:
python plot_genotype_phenotype.py --phenotype Hyperlipidemia \
                                  --binary \
                                  --annotations annotation_points.csv \
                                  --tr-vcf chr1_imputed.vcf.gz
"""

import os
import argparse
import matplotlib.pyplot as plt
import re
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats.contingency import odds_ratio
from gwas_plotter import plot_histogram
from gwas_utils import SAMPLEFILE, GetPTCovarPath, GetCohort, CheckRegion, NormalizeData, \
                        read_annotations, set_genotypes, load_snp_gwas_df, load_gwas_tab

# Scale the size of all texts in the plots including xlabel, ylabel, annotations and legend.
#scale_factor = 2
#plt.rcParams.update({'font.size': plt.rcParams['font.size'] * scale_factor})


def line_plot(data_dict, genotype, phenotype_label, out, confidence_interval=None, std_error=None):
    x = list(data_dict.keys())
    y = list(data_dict.values())
    ax = sns.lineplot(x=x,
                 y=y,
                linewidth=3)
    if confidence_interval:
        ci_x = list(confidence_interval.keys())
        assert len(x) == len(ci_x), f"x: {x} ci_x: {ci_x}"

        lower = [item[0] for item in confidence_interval.values()]
        higher = [item[1] for item in confidence_interval.values()]
        sns.lineplot(x=ci_x, y=lower, alpha=0.5, color='tab:blue')
        sns.lineplot(x=ci_x, y=higher, alpha=0.5, color='tab:blue')
        plt.fill_between(ci_x, lower, higher, alpha=0.5)
    if std_error:
        serr_x = std_error.keys()
        se = std_error.values()
        assert len(x) == len(serr_x), f"len(x): {len(x)} len(serr_x){len(serr_x)}"
        plt.errorbar(x=x, y=y, yerr=list(std_error.values()), fmt='none', c= 'tab:blue')
            
    plt.xlabel(genotype)
    plt.ylabel(phenotype_label)
    plt.savefig(out, bbox_inches="tight")
    plt.clf()

def plot_genotype_phenotype_binary(data, genotype, phenotype, phenotype_label, outpath, epsilon, verbose=False):
    print("Number of cases {} number of controls {} total {}".format(
            len(data[data[phenotype] == 1]),
            len(data[data[phenotype] == 0]),
            len(data),
        ))

    # Check the phenotype is binary
    unique_phenotypes = data[phenotype].unique()
    if len(unique_phenotypes) != 2 or \
            1 not in unique_phenotypes.tolist() or \
            0 not in unique_phenotypes.tolist():
        print("Error: Cannot select --binary if the phenotype is not binary. " + \
              "Unique values for the phenotype {}".format(unique_phenotypes.tolist()))
    
    # Initialize separate variables for each variation of the odds plot    
    data["odds_ratio_threshold"] = None
    odds_ratio_threshold = {}
    odds_ratio_nominal = {}
    odds_threshold = {}
    odds_ratio_threshold_lib = {}
    odds_ratio_threshold_lib_ci = {}
    odds_ratio_threshold_se = {}
    fraction_threshold = {}
    uniq_alleles = sorted(data[genotype].unique())
    min_allele = min(uniq_alleles)

    for uniq_allele in uniq_alleles:
        num_cases_allele = len(data[(data[genotype] == uniq_allele) & (data[phenotype] == 1)]) + epsilon
        num_controls_allele = len(data[(data[genotype] == uniq_allele) & (data[phenotype] == 0)]) + epsilon
        num_cases_no_allele = len(data[(data[genotype] != uniq_allele) & (data[phenotype] == 1)]) + epsilon
        num_controls_no_allele = len(data[(data[genotype] != uniq_allele) & (data[phenotype] == 0)]) + epsilon
        num_cases_gte_allele = len(data[(data[genotype] >= uniq_allele) & (data[phenotype] == 1)]) + epsilon
        num_controls_gte_allele = len(data[(data[genotype] >= uniq_allele) & (data[phenotype] == 0)]) + epsilon
        num_cases_lt_allele = len(data[(data[genotype] < uniq_allele) & (data[phenotype] == 1)]) + epsilon
        num_controls_lt_allele = len(data[(data[genotype] < uniq_allele) & (data[phenotype] == 0)]) + epsilon

        if uniq_allele > min_allele:
            # For threhsold based methods, odds ratio of min value is ill defined.
            odds_ratio_threshold[uniq_allele] = (num_cases_gte_allele/num_cases_lt_allele) \
                                                / (num_controls_gte_allele/num_controls_lt_allele)
            se = np.sqrt(1/num_cases_gte_allele + \
                     1/num_cases_lt_allele + \
                     1/num_controls_gte_allele + \
                     1/num_controls_lt_allele)
            odds_ratio_threshold_se[uniq_allele] = se
            odds_ratio_res = odds_ratio([
                                        [int(num_cases_gte_allele), int(num_cases_lt_allele)],
                                        [int(num_controls_gte_allele), int(num_controls_lt_allele)]])
            odds_ratio_threshold_lib[uniq_allele] = odds_ratio_res.statistic
            odds_ratio_threshold_lib_ci[uniq_allele] = odds_ratio_res.confidence_interval(confidence_level=0.95)
            
            odds_threshold[uniq_allele] = num_cases_gte_allele/num_controls_gte_allele
            fraction_threshold[uniq_allele] = num_cases_gte_allele/(num_controls_gte_allele + num_cases_gte_allele)

        # For nominal method, any allele value can be computed
        odds_ratio_nominal[uniq_allele] = (num_cases_allele/num_cases_no_allele) \
                                            / (num_controls_allele/num_controls_no_allele)
        if verbose:
            print("for gene {} allele {} num_cases(>=allele): {} num_controls(>=allele): {} fraction {}".format(
                  genotype,
                  uniq_allele,
                  num_cases_gte_allele,
                  num_controls_gte_allele,
                  odds_threshold[uniq_allele]))
    
    # Odds ratio nominal plot
    line_plot(data_dict=odds_ratio_nominal,
              genotype=genotype,
              phenotype_label=phenotype_label + " odds ratio == allele",
              out=outpath.replace("genotype", "odds_ratio_nominal"))

    # Odds ratio plot
    line_plot(data_dict=odds_ratio_threshold,
              std_error=odds_ratio_threshold_se,
              genotype=genotype,
              phenotype_label=phenotype_label + " odds ratio >= allele",
              out=outpath.replace("genotype", "odds_ratio_se_error"))

    # Odds ratio plot with error bars
    line_plot(data_dict=odds_ratio_threshold_lib,
              confidence_interval=odds_ratio_threshold_lib_ci,
              genotype=genotype,
              phenotype_label=phenotype_label + " odds ratio >= allele",
              out=outpath.replace("genotype", "odds_ratio_ci_error"))

    # Odds plot
    line_plot(data_dict=odds_threshold,
              genotype=genotype,
              phenotype_label=phenotype_label + " odds >= allele",
              out=outpath.replace("genotype", "odds"))
    
    # Fraction plot
    line_plot(data_dict=fraction_threshold,
              genotype=genotype,
              phenotype_label=phenotype_label + " case control fraction >= allele",
              out=outpath.replace("genotype", "fraction"))

def plot_genotype_phenotype_continuous(data, genotype, phenotype, phenotype_label, outpath, violin, verbose=False):
    # Create a joint grid
    plot = sns.JointGrid()
    # Plot phenotype
    # Plot joint plot in the middle
    if violin:
        sns.violinplot(
            data=data,
            x=genotype,
            y=phenotype,
            ax=plot.ax_joint,
            orient="v",
            native_scale=True,
            )
    else:
        sns.boxplot(
            data=data,
            x=genotype,
            y=phenotype,
            ax=plot.ax_joint,
            orient="v",
            native_scale=True,
            showfliers=False,
            )
    # Plot marginals
    sns.histplot(y=data[phenotype],
                ax=plot.ax_marg_y)
    # Plot genotypes
    sns.histplot(x=data[genotype],
                ax=plot.ax_marg_x)

    plot.set_axis_labels(ylabel=phenotype_label, xlabel=genotype)
    plt.savefig(outpath, bbox_inches="tight")
    plt.clf()

def plot_alleles_histogram(data, genotype, out, epsilon, bin_width=0.5):
    counts = data[genotype].value_counts()
    plot = sns.histplot(data, x=genotype, stat="percent", binwidth=bin_width)
    plot.set_xlabel(genotype)
    # Center the x ticks
    if len(counts) < 10:
        # Center the bars only if there are not too many labels. Otherwise it gets too crowded
        min_val = data[genotype].min()
        max_val = data[genotype].max()
        xtick_labels = np.arange(min_val, max_val + epsilon, bin_width)
        xticks = np.arange(min_val+bin_width/2, max_val+bin_width/2 + epsilon, bin_width)
        plt.xticks(ticks = xticks, labels = xtick_labels)

    plt.savefig(out, bbox_inches="tight")
    plt.clf()


def plot_genotype_phenotype(data, genotype, phenotype,
                            phenotype_label, outpath, binary, violin,
                            merge_bins,
                            min_gt_count=30,
                            min_gt_freq=0.05,
                            verbose=False,
                            epsilon=1E-8):
    plotted_data = data[[genotype, phenotype]].dropna()
    plotted_data = plotted_data.astype(float)

    # Combine bins if indicated to get a more coarse set of alleles by rounding genotypes to the int values
    if merge_bins:
        plotted_data[genotype] = plotted_data[genotype].astype(int)


    # Filter out alleles that are rare
    counts = plotted_data[genotype].value_counts()
    if verbose:
        print("counts before filtering: ", counts)

    # Filter based on counts
    plotted_data = plotted_data.loc[plotted_data[genotype].isin(counts[counts > min_gt_count].index), :]
    if verbose:
        print("counts after count filtering: ", plotted_data[genotype].value_counts())

    # Filter based on frequency
    value_counts = plotted_data[genotype].value_counts()
    to_remove = value_counts[value_counts <= min_gt_freq].index
    plotted_data[genotype] = plotted_data[genotype].replace(to_remove, np.nan)

    # Update the counts after filtering
    counts = plotted_data[genotype].value_counts()
    if verbose:
        print("counts after freq filtering: ", counts)

    # If after filtering we have no polymorphism, return
    if len(counts) == 1:
        print("{} is non-polymorphic after filtering for rare genotypes".format(phenotype))
        return
        
    # Plot alleles histogram after filtering
    out = outpath.replace("genotype", "alleles_hist_after_filter")
    plot_alleles_histogram(data=plotted_data,
                            genotype=genotype,
                            out=out,
                            epsilon=epsilon)

    if binary:
        # Plot a variation of an odds or odds ratio plot.
        plot_genotype_phenotype_binary(data=plotted_data,
                                       genotype=genotype,
                                       phenotype=phenotype,
                                       phenotype_label=phenotype_label,
                                       outpath=outpath,
                                       epsilon=epsilon,
                                       verbose=verbose)
    else:
        # Plot a box plot in the middle and marginal histograms for genotype and phenotype separately.
        plot_genotype_phenotype_continuous(data=plotted_data,
                                           genotype=genotype,
                                           phenotype=phenotype,
                                           phenotype_label=phenotype_label,
                                           outpath=outpath,
                                           violin=violin,
                                           verbose=verbose)

def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes. Required if running associaTR", type=str)
    parser.add_argument("--annotations", help="CSV file with annotations for plotting" + \
                                              "each line should include the " + \
                                              "chromosome, start coordinate and label(gene).",
                                              type=str)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--binary", help="Plot the genotype-phenotype for binary phenotypes", action="store_true")
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--violin", help="Make the genotype_phenotype plot a violinplot instead of boxplot.",
                            action="store_true")
    parser.add_argument("--merge-bins", help="Merge bins in the genotype_phenotype plot.",
                            action="store_true")
    parser.add_argument("--norm", help="Normalize phenotype either quantile or zscore", type=str)
    parser.add_argument("--norm-by-sex",
                        help="Apply the normalization for each sex separately. Default: False",
                        action="store_true")
    parser.add_argument("--outdir", help="Path to the output directory", type=str, required=True)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    # Eventhough covariates are not needed for plotting genotype-phenotypes, we load phenocovar files as
    # sample ids and phenotype lables are stored there.
    ptcovar_path = "phenocovars/{}_phenocovar.csv".format(args.phenotype)
    if os.path.exists(ptcovar_path):
        pass
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)

    # Create output dir if does not exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    # Set up data frame with phenotype and covars
    data = pd.read_csv(ptcovar_path)
    data["person_id"] = data["person_id"].apply(str)

    # Add normalization. If indicated, normalize for each sex separately.
    if args.norm_by_sex:
        # Separate the data into two smaller dataframes based on sex at birth.
        female_data = data[data['sex_at_birth_Male'] == 0].copy()
        male_data = data[data['sex_at_birth_Male'] == 1].copy()

        # Apply normalization on female and male dataframes separately.
        female_data = NormalizeData(data=female_data, norm=args.norm)
        male_data = NormalizeData(data=male_data, norm=args.norm)

        # Concatenate the female and male dataframes back into one
        # and sort the dataframe by original order.
        data = pd.concat([female_data, male_data])
        data = data.sort_index()
    else:
        # Apply normalization on the entire data.
        if args.norm is not None:
            data = NormalizeData(data=data, norm=args.norm)

    # Collect desired samples
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

    # Extract genotypes from vcf and plot allele histogram
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
    plot_histogram(data["phenotype"], args.phenotype,
                os.path.join(args.outdir,
                    "{}_histogram_after_norm.png".format(args.phenotype)))
    
    data = pd.merge(data, samples)

    # Create a more natural phenotype label
    phenotype_label = args.phenotype.lower().replace("_", " ")
    if phenotype_label == "red blood cell distribution width":
        phenotype_label = "RDW"
    if phenotype_label == "malignant neoplasm of the skin":
        phenotype_label = "skin cancer"

    # Plot genotype phenotype
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
                                    gene, args.phenotype, cohort)),
                            min_gt_count=30,
                            min_gt_freq=0.05,
                            verbose=False)

if __name__ == "__main__":
    main()
