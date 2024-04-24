import sys
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="gwas file name", type=str, required=True)
    parser.add_argument("--pval_col", help="name of p value column", type=str, default="P")
    parser.add_argument("--rsid_col", help="name of rsid column", type=str, default="ID")

    args = parser.parse_args()

    gwas = pd.read_csv(args.infile, sep="\t", index_col=False)
    lead_snp_idx = gwas[args.pval_col].argmin()
    lead_snp_id = gwas.iloc[lead_snp_idx][args.rsid_col]
    print(lead_snp_id)




if __name__ == "__main__":
    main()
