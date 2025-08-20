import numpy as np
import pandas as pd
import sys

chrom = sys.argv[1]
lancestry_code = sys.argv[2] #0 for eur

gnomix = f'gnomix-chr{chrom}.msp'

m = np.loadtxt(gnomix, skiprows=2)

with open(gnomix, 'r') as file:
    second_line = file.readlines()[1]
array = np.array(second_line.strip().split('\t'))
columns = list(array)
df = pd.DataFrame(m, columns=columns)
del m
columns = df['spos'].astype('int')
df = df.T
df = df[6:]
df.columns = columns
df = df==0
df = df.astype('int')
df['person_id'] = df.index
df['person_id'] = df['person_id'].apply(lambda x: x.split(".")[0])
a2 = df.groupby('person_id').sum().reset_index()
a2.to_csv(f"chr{chrom}_local_ancestry_{lancestry_code}.csv", index=False)
