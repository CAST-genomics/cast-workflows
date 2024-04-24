import pandas as pd
import sys

f_phenocovar = sys.argv[1]
f_phen = sys.argv[2]
f_covar = sys.argv[3]

phenocovar = pd.read_csv(f_phenocovar)
phen = phenocovar.iloc[:,0:2]
phen.rename(columns={'person_id':'#IID'}, inplace=True)
covar = phenocovar.iloc[:,[0] + list(range(2, len(phenocovar.columns)))] 
covar.rename(columns={'person_id':'IID'}, inplace=True)

phen.to_csv(f_phen, sep="\t", index=False)
covar.to_csv(f_covar, sep="\t", index=False)


