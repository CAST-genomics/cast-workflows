import pandas as pd
import numpy as np
import sys
import scipy.stats as stats

f_samples = sys.argv[1]
f_phen = sys.argv[2]
f_pcs = sys.argv[3]
snppos = sys.argv[4]

def Inverse_Quantile_Normalization(M):
    print("***")
    print(M.shape)
    M = M.transpose()
    R = stats.mstats.rankdata(M)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[0]+1))
    Q = Q.transpose()
    return Q

samples = pd.read_csv(f_samples)
samples.rename(columns={'person_id':'IID'}, inplace=True)
samples = samples['IID']

phen = pd.read_csv(f_phen)
phen.rename(columns={'person_id':'IID'},inplace=True)
pcs = pd.read_csv(f_pcs, sep="\t")
pcs.drop(columns='FID', inplace=True)

snp = np.load(snppos+".npy")
snp = snp.reshape(snp.shape[1],)
geno = pd.DataFrame(data=snp, columns=['snpgeno'])
fam = pd.read_csv(snppos+".fam", sep="\t", header=None)
geno['IID'] = fam[1]
geno = geno[geno['snpgeno'] >= 0]

merge1 = phen.merge(samples, on='IID', how='inner')
merge2 = merge1.merge(pcs, on='IID', how='inner')
merge3 = merge2.merge(geno, on='IID', how='inner')
merge3.insert(0, 'FID', 0)
merge3['phenotype'] = Inverse_Quantile_Normalization(merge3['phenotype'].values)
samples_pre = f_samples[:-4]
merge3.to_csv(f"{samples_pre}_{snppos}", sep="\t", index=False)




