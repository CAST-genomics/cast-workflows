import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def line_plot(data_dict, genotype, phenotype_label, out):
    sns.lineplot(x=data_dict.keys(),
                 y=data_dict.values())
    plt.xlabel(genotype)
    plt.ylabel(phenotype_label + " odds ratio >= allele")
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
    odds_threshold = {}
    fraction_threshold = {}
    uniq_alleles = sorted(data[genotype].unique())

    for uniq_allele in uniq_alleles:
        num_cases_gte_allele = len(data[(data[genotype] >= uniq_allele) & (data[phenotype] == 1)]) + epsilon
        num_controls_gte_allele = len(data[(data[genotype] >= uniq_allele) & (data[phenotype] == 0)]) + epsilon
        num_cases_lt_allele = len(data[(data[genotype] < uniq_allele) & (data[phenotype] == 1)]) + epsilon
        num_controls_lt_allele = len(data[(data[genotype] < uniq_allele) & (data[phenotype] == 0)]) + epsilon
        odds_ratio_threshold[uniq_allele] = (num_cases_gte_allele/num_cases_lt_allele) \
                                            / (num_controls_gte_allele/num_controls_lt_allele)
        odds_threshold[uniq_allele] = num_cases_gte_allele/num_controls_gte_allele
        fraction_threshold[uniq_allele] = num_cases_gte_allele/(num_controls_gte_allele + num_cases_gte_allele)
        if verbose:
            print("for gene {} allele {} num_cases(>=allele): {} num_controls(>=allele): {} fraction {}".format(
                  genotype,
                  uniq_allele,
                  num_cases_gte_allele,
                  num_controls_gte_allele,
                  odds_threshold[uniq_allele]))
    
    # Odds ratio plot
    line_plot(data_dict=odds_ratio_threshold,
              genotype=genotype,
              phenotype_label=phenotype_label + " odds ratio >= allele",
              out=outpath.replace("genotype", "odds_ratio"))

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


def plot_genotype_phenotype(data, genotype, phenotype, gwas, chrom, pos,
                            phenotype_label, outpath, binary, violin,
                            merge_bins,
                            min_gt_count=30,
                            min_gt_freq=0.05,
                            verbose=False,
                            epsilon=1E-8):
    if len(gwas) == 1:
        if gwas.iloc[0]["chrom"] != chrom or \
                abs(int(gwas.iloc[0]["pos"]) - int(pos)) < 1:
            return
    else: # To avoid warning on single element series
        if len(gwas.loc[(gwas["chrom"].astype(str) == chrom) &\
                        (gwas["pos"].astype(int) > int(pos) - 1) &\
                        (gwas["pos"].astype(int) < int(pos) + 1)]) == 0:
            return
    plotted_data = data[[genotype, phenotype]].dropna()
    plotted_data = plotted_data.astype(float)

    # Combine bins if indicated to get a more coarse set of alleles
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
        print("{} is non-polymorphic after filtering for rare genotypes".format(chrom))
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
