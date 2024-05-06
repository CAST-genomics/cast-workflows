import sys
import ast
import pandas as pd

gwas_file = sys.argv[1]
delim="\t"
if len(sys.argv) == 3:
    delim = sys.argv[2]


gwas = pd.read_csv(gwas_file, sep=delim)
alleles = gwas['alleles'].apply(lambda x: ":".join(ast.literal_eval(x)))
gwas['rsid'] = gwas['locus']+":"+alleles
gwas.to_csv(gwas_file, sep="\t", index=False)

