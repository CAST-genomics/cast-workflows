import pandas as pd
import sys
import ast

gwas = pd.read_csv(sys.argv[1], sep="\t")
outfile = sys.argv[2]

gwas['alleles'] = gwas['alleles'].apply(ast.literal_eval)
gwas[['ref', 'alt']] = pd.DataFrame(gwas['alleles'].tolist(), index=gwas.index)
gwas['fam'] = 0

gwas[['chrom', 'rsid', 'fam', 'pos', 'ref', 'alt']].to_csv(outfile, index=False, header=False, sep="\t")
