import pandas as pd
import numpy as np
import sys
import argparse
import re



def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--s1_signals", help="file name for study 1 signals", type=str, required=True)
    parser.add_argument("--s2_signals", help="file name for study 2 signals", type=str, required=True)
    parser.add_argument("--delim", help="delim for signals", type=str, required=True)
    parser.add_argument("--outfile", help="output file name", type=str, required=True)
    parser.add_argument("--pval_col", help="name of p value column", type=str, default="P")
    parser.add_argument("--pos_col", help="name of pos column", type=str, default="POS")
    parser.add_argument("--chr_col", help="name of chrom column", type=str, default="#CHROM")
    parser.add_argument("--id_col", help="name of rsid/rsid equivalent column", type=str, default="ID")
    parser.add_argument("--combine", help="comma separated list of col names to combine for id", type=str, default=None)
    parser.add_argument("--chr_name_prefix", help="if chr col does not contain a number, include prefix here", type=str, default='')

    args = parser.parse_args()

    signals_eur = args.s1_signals
    signals_afr = args.s2_signals
    delim = args.delim
    outfile = args.outfile

    pval_col = args.pval_col
    pos_col = args.pos_col
    chr_col = args.chr_col
    id_col = args.id_col
    id_cols = None
    if id_col == "combine":
        id_col = "combineID"
        assert(args.combine is not None)
        id_cols = [item for item in args.combine.split(",") if item != ""]




    df1 = pd.read_csv(signals_eur, sep=delim)
    print(df1)
    df2 = pd.read_csv(signals_afr, sep=delim)
    print(df2)
    df = pd.concat([df1, df2], ignore_index=True)

    if args.chr_name_prefix != '':
        df[chr_col] = df[chr_col].apply(lambda c: c.removeprefix(args.chr_name_prefix)).astype('int')

    if id_col == "combineID":
        if "alleles" in id_cols:
            df["alleles_concat"] = df["alleles"].apply(lambda allele_lst: '-'.join(re.findall('[AGTC]+', allele_lst)))
            for i in range(len(id_cols)):
                if id_cols[i] == "alleles":
                    id_cols[i] = "alleles_concat"
        df[id_col] = df[id_cols].astype(str).agg('-'.join, axis=1)

    #deduplicate (bc of concat)
    print(df.shape[0])
    df = df.loc[df.groupby(id_col)[pval_col].idxmin()]
    print("DEDUP")
    print(df.shape[0])


    df.sort_values(by = [chr_col, pos_col], inplace=True)

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
        pos = chr_df.loc[snp][pos_col]
        pval = chr_df.loc[snp][pval_col]
        pos_end = pos + w
        pos_start = pos - w
        print(pos)
        print(pos_end)
        print(pos_start)
        nearby = chr_df[ (chr_df[pos_col] <= pos_end) & (chr_df[pos_col] >= pos_start)]
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
