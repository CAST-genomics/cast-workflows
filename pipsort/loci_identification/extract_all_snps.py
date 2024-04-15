import sys
import pandas as pd

gwas1_file = sys.argv[1]
gwas2_file = sys.argv[2]
outfile = sys.argv[3]

gwas1 = pd.read_csv(gwas1_file, sep="\t")
gwas2 = pd.read_csv(gwas2_file, sep="\t")

rsid1 = gwas1['rsid'].values
rsid2 = gwas2['rsid'].values

union_snps = set(rsid1).union(set(rsid2))

with open(outfile, 'w') as f:
    # Iterate over the set and write each element to the file
    for snp in union_snps:
        f.write(snp + '\n')


