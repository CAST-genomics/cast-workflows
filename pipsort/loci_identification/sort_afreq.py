import pandas as pd
import sys

fgwas = sys.argv[1]
fafreq = sys.argv[2]

gwas = pd.read_csv(fgwas, sep="\t")
freq = pd.read_csv(fafreq, sep="\t")
gwas.rename(columns={'rsid':'ID'}, inplace=True)

reordered = freq.set_index("ID").reindex(gwas["ID"]).reset_index()
reordered[freq.columns].to_csv(fafreq, sep="\t", index=False)
