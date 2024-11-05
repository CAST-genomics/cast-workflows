import pandas as pd
import numpy as np
import sys

f_samples = sys.argv[1]
f_phen = sys.argv[2]
f_pcs = sys.argv[3]
f_lancestry = sys.argv[4]
lancestry_code = int(sys.argv[5])
lancestry_pop = sys.argv[6]


lancestry = pd.read_csv(f_lancestry, sep="\t")
print(np.unique(lancestry['eur'], return_counts=True))
lancestry = lancestry[lancestry[lancestry_pop]==lancestry_code]
samples = pd.read_csv(f_samples)
samples = samples.merge(lancestry['person_id'], on='person_id', how='inner')
samples.rename(columns={'person_id': 'IID'}, inplace=True)
print(samples.shape)



phen = pd.read_csv(f_phen)
phen.rename(columns={'person_id':'IID'},inplace=True)
pcs = pd.read_csv(f_pcs, sep="\t")
pcs.drop(columns='FID', inplace=True)

samples = samples['IID']

merge1 = phen.merge(samples, on='IID', how='inner')
merge2 = merge1.merge(pcs, on='IID', how='inner')
merge2.insert(0, 'FID', 0)
samples_pre = f_samples.split(".")[0]
merge2.to_csv(f"{samples_pre}_{lancestry_pop}_{lancestry_code}", sep="\t", index=False)




