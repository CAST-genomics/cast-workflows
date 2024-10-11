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

def plot_histogram(data, outpath):
    plt.clf()
    plot = sns.histplot(data, binwidth=1)
    plt.savefig(outpath)

def plot_genotype_phenotype(data, genotype, phenotype, gwas, chrom, pos, outpath):
    plt.clf()
    plotted_data = data[[genotype, phenotype]].dropna()
    plotted_data = plotted_data.astype(float)
    if len(gwas) == 1:
        if gwas.iloc[0]["chrom"] != chrom or \
                int(gwas.iloc[0]["pos"]) != int(pos):
            print("No gwas data for chrom {} position {}".format(
                chrom, pos))
            return
    else: # To avoid warning on single element series
        if len(gwas.loc[(gwas["chrom"] == chrom) &\
                        (gwas["pos"] == str(pos)), 'beta']) == 0:
            print("No gwas data for chrom {} position {}".format(
                chrom, pos))
            return
    if len(gwas) == 1:
        effect_size = gwas.iloc[0]["beta"]
    else: # To avoid warning on single element series
        effect_size = gwas.loc[(gwas["chrom"] == chrom) &\
                        (gwas["pos"] == str(pos)), 'beta'].item()
    plot = sns.jointplot(
            data=plotted_data,
            alpha=0.5,
            x=genotype,
            y=phenotype)
    print("effect size {:.4}".format(effect_size))
    # Draw a line corresponding to the effect size
    x_0 = data[genotype].mean()
    y_0 = data[phenotype].mean()
    x_1 = data[genotype].mean() - 2 * data[genotype].std()
    x_2 = data[genotype].mean() + 2 * data[genotype].std()

    point_1 = [x_1,
               y_0 + (x_1 - x_0) * effect_size]
    point_2 = [x_2,
               y_0 + (x_2 - x_0) * effect_size]
    plot.ax_joint.plot([point_1[0], x_0, point_2[0]],
                       [point_1[1], y_0, point_2[1]],
                       'b-', linewidth = 1)
    plot.ax_joint.text(x=x_0 + data[genotype].std(),
            y=y_0 + data[phenotype].std(),
            s="effect size: {:.4}".format(effect_size),
            fontsize="medium", weight='bold')
    plt.savefig(outpath)

def annotate_points(ax, gwas):
    # Annotate points
    for point in [("chr15", "88855424", "ACAN"),
                  #("chr17", "30237128", "SLC6A4"),
                  #("chr6", "81752005", "TENT5A"),
                  #("chr20", "2652732", "NOP56"),
                  #("chrX", "43654436", "MAOA"),
                  #("chr15", "101334170", "PCSK6"),
                  #("chr12", "2255790", "CACNA1C"),
                  #("chr1", "155190864", "MUC1"),
                  #("chr21", "43776443", "CSTB"),
                  ]:
        chrom, position, label = point
        locus = gwas[(gwas["chrom"] == chrom) & \
                 (gwas["pos"] == int(position))]
        if len(locus) > 0:
            x = locus["pos"].values[0]
            y = locus["-log10pvalue"].values[0]
            ax.text(x=x, y=y, s=label, fontsize="medium")

def PlotManhattan(gwas, outpath,
                hue,
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
        annotate_points(plot.ax, gwas)

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
