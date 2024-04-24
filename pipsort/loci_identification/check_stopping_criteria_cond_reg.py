import sys
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="gwas file name", type=str, required=True)
    parser.add_argument("--pvalthresh", help="pvalue thresh to stop cond reg", type=float, default=0.001)
    parser.add_argument("--pval_col", help="name of p value column", type=str, default="P")

    args = parser.parse_args()
    
    gwas = pd.read_csv(args.infile, sep="\t", index_col=False)
    stop = int((gwas[args.pval_col] < args.pvalthresh).any())
    print(stop) #0 means stop, 1 means continue


if __name__ == "__main__":
    main()
