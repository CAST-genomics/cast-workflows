import pandas as pd
import numpy as np
import sys
import os


infile = sys.argv[1]
outfile = sys.argv[2]
pval_cutoff = -np.log10(float(sys.argv[3])) #e.g. 5e-08
pval_col = "P"
if len(sys.argv) == 5:
    pval_col = sys.argv[4]

df = pd.read_csv(infile, sep="\t")

df["logp"] = -np.log10(df[pval_col])

df.sort_values(by="logp", inplace=True, ascending=False)

peak_centers = [df.index[0]]

w = 250000
iter = 0
for i in df.index[1:]:
  iter += 1
  if iter % 10000 == 0:
    print("10000 done")
  df_sig = df.loc[peak_centers]
  pos = df.loc[i]['POS']
  pval = df.loc[i]["P"]
  logp = df.loc[i]["logp"]
  if (logp <= pval_cutoff):
    break
  nearby_sig = df_sig[ (df_sig['POS'] <= pos+w) & (df_sig['POS'] >= pos-w)]
  if nearby_sig.shape[0] == 0:
    peak_centers.append(i)

df.loc[peak_centers].sort_values(by=['#CHROM','POS']).to_csv(outfile, index=False)
