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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch
from adjustText import adjust_text



# Scale the size of all texts in the plots including xlabel, ylabel, annotations and legend.
scale_factor = 1.5
plt.rcParams.update({'font.size': plt.rcParams['font.size'] * scale_factor})


def plot_histogram(data, xlabel, outpath, stat='percent', binwidth=0.5):
    plot = sns.histplot(data, binwidth=binwidth, stat=stat)
    plot.set_xlabel(xlabel)
    plt.savefig(outpath, bbox_inches='tight')
    plt.clf()


def annotate_points(ax, gwas, annotations):
    # Annotate points (target TRs) in the manhattan plot based on coordinates.
    texts = []
    for point in annotations:
        chrom, position, label, variant = point
        locus = gwas[(gwas["chrom"] == chrom) & \
                 (gwas["pos"] < int(position) + 1) &\
                 (gwas["pos"] > int(position) - 1)]
        if variant == "snp":
            fontweight = "normal"
        else:
            fontweight = "bold"
        if len(locus) > 1:
            print("multiple VNTRs or associatios within 1 bp for {}:{}. picking the most significant one".format(
                    chrom, position))
        locus = locus.nlargest(1, "-log10pvalue")
        x = locus["pos"].values[0]
        y = locus["-log10pvalue"].values[0]
        texts.append(ax.text(x=x, y=y, s=label, fontsize="medium", fontweight=fontweight))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color='black'))

def PlotManhattan(gwas, outpath,
                hue,
                phenotype,
                annotations,
                annotate=False,
                p_value_threshold=-np.log10(5*10**-8),
                inset_loc="upper center",
                inset_height="35%",
                inset_width="35%",
                inset_window_size_bp=500*1000,
                ):
    num_points = gwas.shape[0]
    size = 30
    chrom = gwas["chrom"][0]
    plotted_gwas = gwas.rename(columns={"-log10pvalue":"-log10pvalue ({})".format(phenotype)})
    plot = sns.relplot(data=gwas, x="pos", y="-log10pvalue",
                s=30, aspect=4, linewidth=0, hue=hue, palette="tab10",
                )
    plot.set_ylabels("-log10pvalue({})".format(phenotype))
    chrom_df = gwas.groupby("chrom")["pos"].median()

    # Set labels
    plot.ax.set_xlabel(chrom + " coordinate")
    # For manhattan plot with more than one chromosome
    #plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.tick_params(axis="both", which="major", labelsize=17)

    # Put the threshold
    plot.ax.axhline(p_value_threshold, linestyle="--", linewidth=1)

    # Set the zommed-in inset
    # Find the window for inset
    max_nlogp = gwas["-log10pvalue"].max()
    tallest_tower_pos = gwas[gwas["-log10pvalue"] == max_nlogp]["pos"].item()
    min_inset_pos = max(tallest_tower_pos - inset_window_size_bp, gwas["pos"].min())
    max_inset_pos = min(tallest_tower_pos + inset_window_size_bp, gwas["pos"].max())
    print("min_inset_pos", min_inset_pos)
    print("max_inset_pos", max_inset_pos)

    # Create an inset dataframe
    inset_gwas = gwas[(gwas["pos"]>=min_inset_pos) & \
                      (gwas["pos"]<=max_inset_pos) & \
                      (gwas["-log10pvalue"] >= p_value_threshold)]
    num_largest = 10
    top_k = inset_gwas.nlargest(num_largest, "-log10pvalue")
    print(f"Top {num_largest} SNP positions based on most significant gwas p-values: ",  top_k)

    # Plot the inset
    colors = {"SNP": "tab:blue", "VNTR": "tab:orange"}
    ax = plot.ax

    axins = inset_axes(ax, width=inset_width, height=inset_height, loc=inset_loc)
    axins.scatter(inset_gwas['pos'], inset_gwas["-log10pvalue"],
                    alpha=0.7,
                    s=20,
                    c=inset_gwas[hue].map(colors))

    axins.set_xlim(min_inset_pos, max_inset_pos)
    axins.set_ylim(p_value_threshold, gwas["-log10pvalue"].max() + 2)
    axins.tick_params(axis="both", which="major", labelsize=17)

    # Remove top and right spines for a more clear look
    axins.spines['top'].set_visible(False)
    axins.spines['right'].set_visible(False)

    # Add a zoomed-in line
    x1, x2 = min_inset_pos, max_inset_pos
    y1, y2 = p_value_threshold, gwas["-log10pvalue"].max()
    
    # In the inset axis
    xB_left, xB_right = axins.get_xlim()
    yB_bottom, yB_top = axins.get_ylim()
    xy_inset1 = (xB_right, yB_bottom)
    xy_inset2 = (xB_right, yB_top)

    # In the main axis
    xy_main1 = (x1, y1)
    xy_main2 = (x2, y2)
    
    # Create the connection patches
    con1 = ConnectionPatch(xyA=xy_inset1, coordsA=axins.transData,
                           xyB=xy_main1, coordsB=ax.transData, color='gray', linestyle='--')
    con2 = ConnectionPatch(xyA=xy_inset2, coordsA=axins.transData,
                           xyB=xy_main2, coordsB=ax.transData, color='gray', linestyle='--')
    # Add to axes
    axins.add_artist(con1)
    axins.add_artist(con2)

    # Annotate selected VNTRs
    if annotate:
        annotate_points(plot.ax, gwas, annotations)
        annotate_points(axins, gwas, annotations)

    plot.fig.savefig(outpath, bbox_inches='tight')

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

        fig.savefig(outpath, bbox_inches='tight')
