import pandas as pd
import argparse
import sys







def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--gwas1", help="gwas s1 file name", type=str, required=True)
    parser.add_argument("--gwas2", help="gwas s2 file name", type=str, required=True)
    parser.add_argument("--fsl", help="sig level for filtering", type=float, required=True)
    parser.add_argument("--pval_col", help="name of pval column", type=str, default="P")
    parser.add_argument("--rsid_col", help="name of rsid column", type=str, default="ID")

    args = parser.parse_args()

    pval_col = args.pval_col
    rsid_col = args.rsid_col
    fsl = args.fsl


    df_s1 = pd.read_csv(args.gwas1, sep="\t")
    df_s2 = pd.read_csv(args.gwas2, sep="\t")

    common = df_s1.merge(df_s2, on=rsid_col, how="inner")
    px = pval_col + "_x"
    py = pval_col + "_y"
    common = common[((common[px] != pd.NA) & (common[px] < fsl)) | ((common[py] != pd.NA) & (common[py] < fsl))]
    common = common[rsid_col]

    df_s1_sig = df_s1.merge(common, on=rsid_col, how="inner")

    df_s2_sig = df_s2.merge(common, on=rsid_col, how="inner")

    assert(df_s1_sig.shape[0] == df_s2_sig.shape[0])

    df_s1_sig.to_csv("s1_gwas_temp.txt", sep="\t", index=False)
    df_s2_sig.to_csv("s2_gwas_temp.txt", sep="\t", index=False)

    df_s1_sig[rsid_col].to_csv("s1_snps.txt", header=False, index=False)
    df_s2_sig[rsid_col].to_csv("s2_snps.txt", header=False, index=False)


if __name__ == "__main__":
    main()
