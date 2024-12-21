import pandas as pd
import sys

fregions = sys.argv[1]
outfile = sys.argv[2]
pop = sys.argv[3]

df = pd.read_csv(fregions, sep="\t")
df2 = df.T
df3 = df2[6:]
df3[0] = (df3[0]==0).astype('int')
df3['person_id'] = df3.index
df3['person_id'] = df3['person_id'].apply(lambda x: x.split(".")[0])
a2 = df3.groupby('person_id')[0].sum().reset_index()
a2.rename(columns={0:pop}, inplace=True)
a2.to_csv(f"{outfile}", sep="\t", index=False)
