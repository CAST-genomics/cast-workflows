import sys
import ast
import pandas as pd

gwas_file = sys.argv[1]

gwas = pd.read_csv(gwas_file, sep="\t")
alleles = gwas['alleles'].apply(lambda x: ":".join(ast.literal_eval(x)))
gwas['rsid'] = gwas['locus']+":"+alleles
gwas.to_csv(gwas_file, sep="\t", index=False)

