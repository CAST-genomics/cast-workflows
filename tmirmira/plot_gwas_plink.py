#from cast-workflows/gwas/aou
import matplotlib.pyplot as plt
import numbers
import numpy as np
import pandas as pd
import seaborn as sns




def PlotManhattan(gwasfile, outpath):
    gwas = pd.read_csv(gwasfile, sep="\t")
    gwas["-log10pvalue"] = np.where(
        gwas['P'] == 0,
        300,
        -np.log10(gwas['P'])
    )
    gwas["ind"] = range(gwas.shape[0])
    plot = sns.relplot(data=gwas, x="ind", y="-log10pvalue", \
        s=6, aspect=4, linewidth=0, hue="#CHROM", palette="tab10", legend=None)
    chrom_df = gwas.groupby("#CHROM")["ind"].median()
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

def PlotQQ(gwasfile, outpath):
    gwas = pd.read_csv(gwasfile, sep="\t")
    gwas["-log10pvalue"] = np.where(
        gwas['P'] == 0,
        300,
        -np.log10(gwas['P'])
    )
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

if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    outputfilename = sys.argv[2]
    PlotManhattan(filename, f'gwas_{outputfilename}')
    PlotQQ(filename, f'qq_{outputfilename}')
