import pandas as pd
import numpy as np
import sys
import scipy.stats as stats

f_samples = sys.argv[1]
f_phen = sys.argv[2]
f_pcs = sys.argv[3]
f_lancestry = sys.argv[4]
lancestry_code = int(sys.argv[5])
lancestry_pop = sys.argv[6]

def Inverse_Quantile_Normalization(M):
    print("***")
    print(M.shape)
    M = M.transpose()
    R = stats.mstats.rankdata(M)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[0]+1))
    Q = Q.transpose()
    return Q

lancestry = pd.read_csv(f_lancestry, sep="\t")
print(np.unique(lancestry[lancestry_pop], return_counts=True))
if lancestry_code >= 0: # filter. If negative, then don't filter
    lancestry = lancestry[lancestry[lancestry_pop]==lancestry_code]
samples = pd.read_csv(f_samples)
if lancestry_code >= 0:
    samples = samples.merge(lancestry['person_id'], on='person_id', how='inner')
else:
    samples = samples.merge(lancestry, on='person_id', how='inner')
samples.rename(columns={'person_id': 'IID'}, inplace=True)
print(samples.shape)



phen = pd.read_csv(f_phen)
phen.rename(columns={'person_id':'IID'},inplace=True)
pcs = pd.read_csv(f_pcs, sep="\t")
pcs.drop(columns='FID', inplace=True)

if lancestry_code >= 0:
    samples = samples['IID']
else:
    samples = samples[['IID', lancestry_pop]]

merge1 = phen.merge(samples, on='IID', how='inner')
merge2 = merge1.merge(pcs, on='IID', how='inner')
merge2.insert(0, 'FID', 0)
#merge2['phenotype'] = Inverse_Quantile_Normalization(merge2['phenotype'].values)
#samples_pre = f_samples.split(".")[0]
samples_pre = f_samples[:-4]
if lancestry_code >= 0:
    merge2.to_csv(f"{samples_pre}_{lancestry_pop}_{lancestry_code}", sep="\t", index=False)
else: #we didn't filter on lancestry code so write out diff file name
    merge2.to_csv(f"{samples_pre}_{lancestry_pop}", sep="\t", index=False)




