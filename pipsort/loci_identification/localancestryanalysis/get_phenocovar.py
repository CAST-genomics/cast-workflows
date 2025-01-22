import pandas as pd
import numpy as np
import sys
import scipy.stats as stats

f_samples = sys.argv[1]
f_phen = sys.argv[2]
f_pcs = sys.argv[3]
outfile = sys.argv[4]

def Inverse_Quantile_Normalization(M):
    print("***")
    print(M.shape)
    M = M.transpose()
    R = stats.mstats.rankdata(M)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[0]+1))
    Q = Q.transpose()
    return Q

samples = pd.read_csv(f_samples)
samples.rename(columns={'research_id': 'IID'}, inplace=True)


phen = pd.read_csv(f_phen)
phen.rename(columns={'person_id':'IID'},inplace=True)
pcs = pd.read_csv(f_pcs, sep="\t")
pcs.drop(columns='FID', inplace=True)


merge1 = phen.merge(samples, on='IID', how='inner')
merge2 = merge1.merge(pcs, on='IID', how='inner')
merge2.insert(0, 'FID', 0)
merge2['phenotype'] = Inverse_Quantile_Normalization(merge2['phenotype'].values)
merge2.to_csv(outfile, sep="\t", index=False)
