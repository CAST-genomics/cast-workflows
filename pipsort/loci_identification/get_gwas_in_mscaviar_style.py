import argparse
import pandas as pd
import numpy as np
import sys


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="gwas file name", type=str, required=True)
    parser.add_argument("--outfile", help="output file name", type=str, required=True)
    parser.add_argument("--chr_col", help="name of chr column", type=str, default="#CHROM")
    parser.add_argument("--pos_col", help="name of pos column", type=str, default="POS")
    parser.add_argument("--rsid_col", help="name of rsid column", type=str, default="ID")
    parser.add_argument("--ref_col", help="name of ref column", type=str, default="REF")
    parser.add_argument("--alt_col", help="name of alt column", type=str, default="ALT")
    parser.add_argument("--beta_col", help="name of beta column", type=str, default="BETA")
    parser.add_argument("--se_col", help="name of se column", type=str, default="SE")

    args = parser.parse_args()

    chr_col = args.chr_col
    pos_col = args.pos_col
    rsid_col = args.rsid_col
    ref_col = args.ref_col
    alt_col = args.alt_col
    beta_col = args.beta_col
    se_col = args.se_col


    gwas = pd.read_csv(args.infile, sep="\t", comment='#')
    gwas = gwas[[chr_col, pos_col, ref_col, alt_col, rsid_col, beta_col, se_col]]
    gwas.rename(columns={chr_col:'CHROMOSOME',pos_col:'BP', ref_col: 'REF_ALLELE', alt_col: 'ALT_ALLELE', rsid_col:'SNP_ID', beta_col:'BETA', se_col:'SE' }, inplace=True)

    gwas.to_csv(args.outfile, sep="\t", index=False, na_rep='NA')
    gwas = pd.read_csv(args.outfile, sep="\t")
    gwas.dropna(axis=0, inplace=True)
    gwas.to_csv(args.outfile, sep="\t", index=False, na_rep='NA')

if __name__ == "__main__":
    main()
