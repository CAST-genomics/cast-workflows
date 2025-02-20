"""
Make Manhattan and QQ plots
QQ plot code copied from:
https://github.com/satchellhong/qqman/
"""

import matplotlib.pyplot as plt
import numbers
import numpy as np
import pandas as pd
import seaborn as sns

def plot_histogram(data, xlabel, outpath):
    plot = sns.histplot(data, binwidth=0.5)
    plot.set_xlabel(xlabel)
    plt.savefig(outpath)
    plt.clf()

def plot_genotype_phenotype(data, genotype, phenotype, gwas, chrom, pos,
                            phenotype_label, outpath, downsample, violin,
                            merge_bins,
                            min_gt_count=10):
    if len(gwas) == 1:
        if gwas.iloc[0]["chrom"] != chrom or \
                abs(int(gwas.iloc[0]["pos"]) - int(pos)) < 1:
            return
    else: # To avoid warning on single element series
        if len(gwas.loc[(gwas["chrom"].astype(str) == chrom) &\
                        (gwas["pos"].astype(int) > int(pos) - 1) &\
                        (gwas["pos"].astype(int) < int(pos) + 1)]) == 0:
            return
    if len(gwas) == 1:
        effect_size = gwas.iloc[0]["beta"]
    else: # To avoid warning on single element series
        effect_size = gwas.loc[(gwas["chrom"] == chrom) &\
                        (gwas["pos"] > int(pos) - 1) &\
                        (gwas["pos"] < int(pos) + 1), 'beta'].item()
    plotted_data = data[[genotype, phenotype]].dropna()
    plotted_data = plotted_data.astype(float)

    if downsample:
        # This is only necessary for the binary phenotypes.
        # Where usually negative samples are over represented and 
        # should be sampled
        # Check the phenotype is binary
        unique_phenotypes = plotted_data[phenotype].unique()
        if len(unique_phenotypes) != 2 or \
                1 not in unique_phenotypes.tolist() or \
                0 not in unique_phenotypes.tolist():
            print("Error: Cannot select --downsample if the phenotype is not binary. " + \
                  "Unique values for the phenotype {}".format(unique_phenotypes.tolist()))

            
        num_cases = len(plotted_data[plotted_data[phenotype] == 1])
        num_controls = len(plotted_data[plotted_data[phenotype] == 0])
        cases = plotted_data[plotted_data[phenotype] == 1]
        controls_ds = plotted_data[plotted_data[phenotype] == 0].sample(n=num_cases)
        plotted_data = pd.concat([controls_ds, cases])
        print("Down sampling {} controls to {} in order to match {} cases".format(
                    num_controls, len(controls_ds), num_cases))

    # Create a joint grid
    plot = sns.JointGrid()
    if merge_bins:
        plotted_data[genotype] = plotted_data[genotype].astype(int)
    counts = plotted_data[genotype].value_counts()
    print("counts before filtering: ", counts)
    plotted_data = plotted_data.loc[plotted_data[genotype].isin(counts[counts > min_gt_count].index), :]
    print("counts after filtering: ", plotted_data[genotype].value_counts())
    # Plot phenotype
    # Plot joint plot in the middle
    if violin:
        sns.violinplot(
            data=plotted_data,
            x=genotype,
            y=phenotype,
            ax=plot.ax_joint,
            orient="v",
            native_scale=True,
            )
    else:
        sns.boxplot(
            data=plotted_data,
            x=genotype,
            y=phenotype,
            ax=plot.ax_joint,
            orient="v",
            native_scale=True,
            showfliers=False,
            )
    # Plot marginals
    sns.histplot(y=plotted_data[phenotype],
                ax=plot.ax_marg_y)
    # Plot genotypes
    sns.histplot(x=plotted_data[genotype],
                ax=plot.ax_marg_x)

    plot.set_axis_labels(ylabel=phenotype_label, xlabel=genotype)
    print("effect size {:.4}".format(effect_size))
    plt.savefig(outpath, bbox_inches="tight")
    plt.clf()

def annotate_points(ax, gwas, annotations):
    # Annotate points
    for point in annotations:
        chrom, position, label = point
        locus = gwas[(gwas["chrom"] == chrom) & \
                 (gwas["pos"] < int(position) + 1) &\
                 (gwas["pos"] > int(position) - 1)]
        if len(locus) > 0:
            print("multiple VNTRs within 1 bp for {}:{}".format(chrom, position))
            x = locus["pos"].values[0]
            y = locus["-log10pvalue"].values[0]
            ax.text(x=x, y=y, s=label, fontsize="medium")

def PlotManhattan(gwas, outpath,
                hue,
                annotations,
                annotate=False,
                p_value_threshold=-np.log10(5*10**-8),
                extra_points=None
                ):
    num_points = gwas.shape[0]
    size = 30
    plot = sns.relplot(data=gwas, x="pos", y="-log10pvalue",
                s=30, aspect=4, linewidth=0, hue=hue, palette="tab10",
                )
    chrom_df = gwas.groupby("chrom")["pos"].median()

    if extra_points and num_points > 0 and max(gwas['-log10pvalue'].dropna()):
        for chrom, pos, p_value, label in extra_points:
            chrom_gwas = gwas[gwas["chrom"] == chrom] 
            # Index of the new node in the manhattan plot is the number of values
            # with position smaller than the current position.
            index = len(chrom_gwas[chrom_gwas["pos"] < pos])
            plot.axes[0][0].scatter(x=index, y=p_value, color="red", marker="^")
            plot.axes[0][0].text(x=index, y=p_value, s=label, fontsize="medium")
            plot.set(ylim=(0, max(gwas['-log10pvalue'].dropna()) + 1), xlim=(0, num_points))


    if annotate:
        annotate_points(plot.ax, gwas, annotations)

    # Set labels
    plot.ax.set_xlabel("Chromosome")
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)

    # Put the threshold
    plot.ax.axhline(p_value_threshold, linestyle="--", linewidth=1)
    plot.fig.savefig(outpath)

def ppoints(n, a=None):
        """ numpy analogue or `R`'s `ppoints` function
                see details at https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/ppoints
                https://docs.tibco.com/pub/enterprise-runtime-for-R/5.0.0/doc/html/Language_Reference/stats/ppoints.html
                :param n: array type or number"""
        
        if isinstance(n, numbers.Number):
                n = float(n)
        else:
                n = float(len(n))
        if a == None:
                a = .375 if n<=10 else .5
                
        return (np.arange(n) + 1 - a)/(n + 1 - 2*a)

def PlotQQ(gwas, outpath):
        p_vals = gwas["-log10pvalue"].dropna().apply(lambda x: 10**-x)
        p_vals = p_vals[(0<p_vals)&(p_vals<1)]
        observed = -np.log10(np.sort(np.array(p_vals)))
        expected = -np.log10(ppoints(len(p_vals)))
        x_padding = (np.nanmax(expected)-np.nanmin(expected))/12
        y_padding = (np.nanmax(observed)-np.nanmin(observed))/12

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (12,12))

        ax.scatter(expected,observed,c='k')
        
        xlim_min = np.nanmin(expected)-x_padding
        xlim_max = np.nanmax(expected)+x_padding
        ylim_min = np.nanmin(observed)-y_padding
        ylim_max = np.nanmax(observed)+y_padding
        
        max_lim = xlim_max if xlim_max<ylim_max else ylim_max
        min_lim = xlim_min if xlim_min>ylim_min else ylim_min
        ax.plot([min_lim,max_lim],[min_lim,max_lim],'r-')

        if xlim_min < xlim_max and ylim_min < ylim_max:
            ax.set_xlim([xlim_min, xlim_max])
            ax.set_ylim([ylim_min, ylim_max])
        ax.set_xlabel("Expected $-log_{10}(p)$")
        ax.set_ylabel("Observed $-log_{10}(p)$")

        fig.savefig(outpath)
