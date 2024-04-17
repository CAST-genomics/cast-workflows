import argparse
import pandas as pd
import numpy as np
import sys


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="gwas file name", type=str, required=True)
    parser.add_argument("--pos_col", help="name of pos column", type=str, default="POS")
    parser.add_argument("--rsid_col", help="name of rsid column", type=str, default="ID")

    args = parser.parse_args()

    pos_col = args.pos_col
    rsid_col = args.rsid_col


    gwas = pd.read_csv(args.infile, sep="\t")

    gwas.sort_values(pos_col, inplace=True, ignore_index=True)

    orig_idx = np.arange(gwas.shape[0], dtype=int)
    curr_pos = gwas[pos_col][0]
    curr_idx = 0
    for i in range(1, gwas.shape[0]):
        next_pos = gwas[pos_col][i]
        if (next_pos != curr_pos):
            curr_idx += 1
        orig_idx[i] = curr_idx
        curr_pos = next_pos
    gwas["orig_idx"] = orig_idx
    gwas.sort_values(["orig_idx",rsid_col], inplace=True, ignore_index=True)
    gwas.drop(columns=["orig_idx"],inplace=True)
    gwas.to_csv(args.infile, sep="\t", index=False)

if __name__ == "__main__":
    main()
