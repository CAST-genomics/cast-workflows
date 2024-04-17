import pandas as pd
import argparse
import sys







def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--gwas1", help="gwas s1 file name", type=str, required=True)
    parser.add_argument("--gwas2", help="gwas s2 file name", type=str, required=True)
    parser.add_argument("--fsl1", help="sig level for filtering s1", type=float, required=True)
    parser.add_argument("--fsl2", help="sig level for filtering s2", type=float, default=None)
    parser.add_argument("--pval_col", help="name of pval column", type=str, default="P")
    parser.add_argument("--rsid_col", help="name of rsid column", type=str, default="ID")

    args = parser.parse_args()

    pval_col = args.pval_col
    rsid_col = args.rsid_col
    if args.fsl2 is None:
        args.fsl2 = args.fsl1
    fsl1 = args.fsl1
    fsl2 = args.fsl2


    df_s1 = pd.read_csv(args.gwas1, sep="\t")
    df_s2 = pd.read_csv(args.gwas2, sep="\t")

    df_s1_sig = df_s1.merge(df_s2, on=rsid_col, how="left")
    p1 = pval_col+"_x"
    p2 = pval_col+"_y"
    df_s1_sig = df_s1_sig[((df_s1_sig[p1] != pd.NA) & (df_s1_sig[p1] < fsl1)) | ((df_s1_sig[p2] != pd.NA) & (df_s1_sig[p2] < fsl2))]

    df_s2_sig = df_s2.merge(df_s1, on=rsid_col, how="left")
    df_s2_sig = df_s2_sig[((df_s2_sig[p1] != pd.NA) & (df_s2_sig[p1] < fsl1)) | ((df_s2_sig[p2] != pd.NA) & (df_s2_sig[p2] < fsl2))]

    df_s1_sig_final = df_s1.merge(df_s1_sig[rsid_col], on=rsid_col, how="inner")
    df_s2_sig_final = df_s2.merge(df_s2_sig[rsid_col], on=rsid_col, how="inner")

    df_s1_sig_final.to_csv("s1_gwas_temp.txt", sep="\t", index=False)
    df_s2_sig_final.to_csv("s2_gwas_temp.txt", sep="\t", index=False)

    df_s1_sig_final[rsid_col].to_csv("s1_snps.txt", header=False, index=False)
    df_s2_sig_final[rsid_col].to_csv("s2_snps.txt", header=False, index=False)

if __name__ == "__main__":
    main()
