import pandas as pd
import sys
import ast

fgwas = sys.argv[1]
gwas = pd.read_csv(fgwas, sep="\t")

gwas['alleles'] = gwas['alleles'].apply(ast.literal_eval)
gwas[['ref', 'alt']] = pd.DataFrame(gwas['alleles'].tolist(), index=gwas.index)

def get_a1(row):
    if row['beta'] >= 0:
        return row['ref']
    else:
        return row['alt']
def get_a2(row):
    if row['beta'] >= 0:
        return row['alt']
    else:
        return row['ref']

gwas['a1'] = gwas.apply(get_a1, axis=1)
gwas['a2'] = gwas.apply(get_a2, axis=1)
chrom = gwas['chrom'].iloc[0][3:]
gwas['chrom'] = chrom


gwas[['chrom', 'rsid', 'pos', 'ref', 'alt', 'beta', 'standard_error', 'p_value']].to_csv(fgwas, index=False, sep="\t")
