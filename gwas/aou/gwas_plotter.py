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



def plot_histogram(data, xlabel, outpath, stat='percent', binwidth=0.5):
    plot = sns.histplot(data, binwidth=binwidth, stat=stat)
    plot.set_xlabel(xlabel)
    plt.savefig(outpath, bbox_inches='tight')
    plt.clf()


def annotate_points(ax, gwas, annotations):
    # Annotate points (target TRs) in the manhattan plot based on coordinates.
    for point in annotations:
        chrom, position, label = point
        locus = gwas[(gwas["chrom"] == chrom) & \
                 (gwas["pos"] < int(position) + 1) &\
                 (gwas["pos"] > int(position) - 1)]
        if len(locus) > 1:
            print("multiple TRs or associatios within 1 bp for {}:{}".format(chrom, position))
        elif len(locus) == 1:
            x = locus["pos"].values[0]
            y = locus["-log10pvalue"].values[0]
            ax.text(x=x, y=y, s=label, fontsize="medium")


def PlotManhattan(gwas, outpath):
    gwas["ind"] = range(gwas.shape[0])
    plot = sns.relplot(data=gwas, x="ind", y="-log10pvalue", \
        s=6, aspect=4, linewidth=0, hue="chrom", palette="tab10", legend=None)
    chrom_df = gwas.groupby("chrom")["ind"].median()
    plot.ax.set_xlabel("Chromosome")
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.axhline(-np.log10(5*10**-8), linestyle="--", linewidth=1)
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
