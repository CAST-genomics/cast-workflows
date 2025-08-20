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
