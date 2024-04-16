import pandas as pd
import sys

phen_file = sys.argv[1]
outfile = sys.argv[2]

phen = pd.read_csv(phen_file)
FID = pd.DataFrame(0, index=phen.index, columns=['FID'])
FID['IID'] = phen['person_id']
FID['PHEN'] = phen['phenotype']
FID.to_csv(outfile, sep="\t", index=False)
