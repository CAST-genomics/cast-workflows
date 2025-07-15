import argparse
import pandas as pd
import numpy as np
import sys
import os


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="gwas file name", type=str, required=True)
    parser.add_argument("--outfile", help="output file name", type=str, required=True)
    parser.add_argument("--pval_cutoff", help="sig pvalue cutoff", type=float, default=5e-08)
    parser.add_argument("--pval_col", help="name of p value column", type=str, default="P")
    parser.add_argument("--pos_col", help="name of pos column", type=str, default="POS")
    parser.add_argument("--chr_col", help="name of chrom column", type=str, default="#CHROM")
    parser.add_argument("--comment", help="comment symbol for header lines", type=str, default=None)

    args = parser.parse_args()

    logpval_cutoff = -np.log10(args.pval_cutoff)
    pval_col = args.pval_col
    pos_col = args.pos_col
    chr_col = args.chr_col
    cmt = args.comment


    df = pd.read_csv(args.infile, sep="\t", comment=cmt)

    df["logp"] = -np.log10(df[pval_col])
    df.sort_values(by="logp", inplace=True, ascending=False) 


    infloc = -1
    if df.iloc[0]['logp'] == np.inf:
        infloc = 0
        while df.iloc[infloc]['logp'] == np.inf:
            infloc += 1
        infdf = df.iloc[:infloc]
        df.iloc[:infloc] = infdf.sort_values(by=[chr_col, pos_col], ascending=True)
    

    peak_centers = [df.index[0]]

    w = 250000
    iter = 0
    for i in df.index[1:]:
        iter += 1
        if iter % 10000 == 0:
            print("10000 done")
        df_sig = df.loc[peak_centers]
        pos = df.loc[i][pos_col]
        pval = df.loc[i][pval_col]
        logp = df.loc[i]["logp"]
        if (logp <= logpval_cutoff):
            break
        nearby_sig = df_sig[ (df_sig[pos_col] <= pos+w) & (df_sig[pos_col] >= pos-w)]
        if nearby_sig.shape[0] == 0:
            peak_centers.append(i)

    df.loc[peak_centers].sort_values(by=[chr_col,pos_col]).to_csv(args.outfile, index=False)

if __name__ == "__main__":
    main()