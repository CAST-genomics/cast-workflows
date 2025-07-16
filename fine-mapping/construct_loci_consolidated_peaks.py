import pandas as pd
import numpy as np
import sys
import argparse
import re



def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--signals", help="file name for all signals", type=str, required=True)
    parser.add_argument("--outfile", help="output file name", type=str, required=True)
    parser.add_argument("--delim", help="delim for signals", type=str, default=",")
    parser.add_argument("--pos_col", help="name of pos column", type=str, default="POS")
    parser.add_argument("--chr_col", help="name of chrom column", type=str, default="#CHROM")
    parser.add_argument("--pval_col", help="name of pval column", type=str, default="P")
    parser.add_argument("--id_col", help="name of rsid/rsid equivalent column", type=str, default="ID")
    parser.add_argument("--chr_name_prefix", help="if chr col does not contain a number, include prefix here", type=str, default='')

    args = parser.parse_args()

    signals = args.signals
    outfile = args.outfile

    delim = args.delim
    pos_col = args.pos_col
    chr_col = args.chr_col
    pval_col = args.pval_col
    id_col = args.id_col

    df = pd.read_csv(signals, sep=delim)

    if args.chr_name_prefix != '':
        df[chr_col] = df[chr_col].apply(lambda c: c.removeprefix(args.chr_name_prefix)).astype('int')


    df.sort_values(by = [chr_col, pos_col], inplace=True)

    #exclude snps in mhc region chr6:28,510,120-33,480,577
    df = df[~((df[chr_col]==6)&(df[pos_col]>=28510120)&(df[pos_col]<=33480577))]

    #deduplicate (bc of concat)
    print(df.shape[0])
    df = df.loc[df.groupby(id_col)[pval_col].idxmin()]
    print("DEDUP")
    print(df.shape[0])

    loci = []
    w = 1000000
    w2 = w / 2
    for i in range(1,23): #22 chromosomes
      print(i)
      chr_df = df[df[chr_col] == i]
      snps = list(chr_df[id_col].values)
      print(snps)
      chr_df.set_index(id_col, inplace=True)
      print(chr_df)
      while len(snps) > 0:
        num_snps = 1
        snp = snps.pop(0)
        print("snp ", snp)
        pos = chr_df.loc[snp][pos_col]
        pos_end = pos + w
        pos_start = pos - w
        print("pos ", pos)
        print("pos end ", pos_end)
        print("pos start ", pos_start)
        nearby = chr_df[ (chr_df[pos_col] <= pos_end) & (chr_df[pos_col] >= pos_start)]
        print("nearby")
        print(nearby)
        region = [snp]
        while nearby.shape[0] > num_snps:
          new_snps = list(set(nearby.index).difference(set(region)))
          for s in new_snps:
            snps.remove(s)
          region += new_snps
          num_snps = len(region)
          pos_end = np.max(nearby[pos_col]+w)
          pos_start = np.min(nearby[pos_col]-w)
          nearby = chr_df[ (chr_df[pos_col] <= pos_end) & (chr_df[pos_col] >= pos_start)]


        loci.append([i, pos_start+w2, pos_end-w2, num_snps])
        print([i, pos_start+w2, pos_end-w2])

    loci_arr = np.array(loci, dtype=int)
    np.savetxt(outfile, loci_arr, fmt="%s")

if __name__ == "__main__":
    main()