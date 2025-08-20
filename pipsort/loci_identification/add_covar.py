import sys
import pandas as pd
import re
import argparse



def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--covars", help="gwas file name", type=str, required=True)
    parser.add_argument("--genofile", help="genotype to add as covar file name", type=str, required=True)
    parser.add_argument("--snp", help="name of snp to add as covar", type=str, required=True)

    args = parser.parse_args()

    covars = pd.read_csv(args.covars, sep="\t", index_col=False)
    genotypes = pd.read_csv(args.genofile, sep="\t", index_col=False)
    snp_name = args.snp

    columns = genotypes.columns
    r = re.compile(snp_name+"*")
    newlist = list(filter(r.match, columns))
    assert(len(newlist) == 1)
    col_name = newlist[0]


    genotypes = genotypes[['IID', col_name]]

    covars = covars.merge(genotypes, how="left", left_on=['IID'], right_on=['IID'])

    covars.to_csv(args.covars, sep="\t", index=False, na_rep="NA")

if __name__ == "__main__":
    main()
