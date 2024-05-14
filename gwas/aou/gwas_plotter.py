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
    plot = sns.histplot(data)
    plt.savefig(outpath)

def plot_genotype_phenotype(data, genotype, phenotype, gwas, chrom, pos, outpath):
    plt.clf()
    plotted_data = data[[genotype, phenotype]].dropna()
    plotted_data = plotted_data.astype(float)
    plot = sns.jointplot(
            data=plotted_data,
            alpha=0.5,
            x=genotype,
            y=phenotype)
    effect_size = gwas.loc[(gwas["chrom"] == chrom) &\
                        (gwas["pos"] == pos), 'beta'].item()
    #print("effect size {:.4}".format(effect_size))
    # Draw a line corresponding to the effect size
    #point_1 = [min(data[genotype]), min(data[phenotype])]
    x_0 = data[genotype].mean()
    y_0 = data[phenotype].mean()
    x_1 = data[genotype].mean() - 2 * data[genotype].std()
    x_2 = data[genotype].mean() + 2 * data[genotype].std()

    point_1 = [x_1,
               #data[phenotype].mean() - data[phenotype].std()]
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
                  ("chr17", "30237128", "SLC6A4"),
                  ("chr6", "81752005", "TENT5A"),
                  ("chr20", "2652732", "NOP56"),
                  ("chrX", "43654436", "MAOA"),
                  ("chr15", "101334170", "PCSK6"),
                  ("chr12", "2255790", "CACNA1C"),
                  ("chr1", "155190864", "MUC1"),
                  ("chr21", "43776443", "CSTB"),
                  ]:
        chrom, position, label = point
        locus = gwas[(gwas["chrom"] == chrom) & \
                 (gwas["pos"] == int(position))]
        x = locus["ind"].values[0]
        y = locus["-log10pvalue"].values[0]
        ax.text(x=x, y=y, s=label, fontsize="medium")

def PlotManhattan(gwas, outpath):
    gwas["ind"] = range(gwas.shape[0])
    plot = sns.relplot(data=gwas, x="ind", y="-log10pvalue", \
        s=30, aspect=4, linewidth=0, hue="chrom", palette="tab10", legend=None)
    chrom_df = gwas.groupby("chrom")["ind"].median()

    annotate_points(plot.ax, gwas)

    # Set labels
    plot.ax.set_xlabel("Chromosome")
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)

    # Put the threshold
    #plot.ax.axhline(-np.log10(5*10**-8), linestyle="--", linewidth=1)
    plot.ax.axhline(-np.log10(1*10**-3), linestyle="--", linewidth=1)
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
        
        ax.set_xlim([xlim_min, xlim_max])
        ax.set_ylim([ylim_min, ylim_max])
        ax.set_xlabel("Expected $-log_{10}(p)$")
        ax.set_ylabel("Observed $-log_{10}(p)$")

        fig.savefig(outpath)
