import pandas as pd
import sys

s1_samples = pd.read_csv(sys.argv[1])
s2_samples = pd.read_csv(sys.argv[2])
phen = pd.read_csv(sys.argv[3])
f_s1gwas = sys.argv[4]
f_s2gwas = sys.argv[5]

s1_samples = s1_samples.merge(phen, on='person_id', how='inner')
s2_samples = s2_samples.merge(phen, on='person_id', how='inner')

s1gwas = pd.read_csv(f_s1gwas, sep="\t")
s1gwas['n'] = s1_samples.shape[0]
s1gwas.to_csv(f_s1gwas, sep="\t", index=False)

s2gwas = pd.read_csv(f_s2gwas, sep="\t")
s2gwas['n'] = s2_samples.shape[0]
s2gwas.to_csv(f_s2gwas, sep="\t", index=False)

