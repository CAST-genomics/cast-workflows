import numpy as np
import pandas as pd
import sys

chrom = sys.argv[1]

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

labels={0:'EUR',1:'EAS' ,2:'NAT',3:'AFR', 4:'SAS', 5:'AHG' ,6:'OCE',7:'WAS'}
labels['index']='research_id'

results = {}
for col in df.columns:
    unique_values = np.unique(df[col], return_counts=True)[1]
    results[col] = unique_values

result_df = pd.DataFrame.from_dict(results, orient='index')
result_df.reset_index(inplace=True)
result_df.rename(columns=labels, inplace=True)
result_df.to_csv(f"chr{chrom}_count_ancestry_haplotypes_per_region.csv", index=False)
row_sum = result_df.iloc[0, 1:].sum()
result_df.iloc[:, 1:] = result_df.iloc[:, 1:].div(row_sum)
result_df.to_csv(f"chr{chrom}_percent_ancestry_haplotypes_per_region.csv", index=False)



df2 = df.T
unique_values_dict = df2.apply(lambda col: pd.Series(np.unique(col, return_counts=True)[1]), axis=0)
result_df = unique_values_dict.T.reset_index()
result_df.rename(columns=labels, inplace=True)
result_df['research_id'] = result_df['research_id'].apply(lambda x: x.split(".")[0])
result_df.fillna(0, inplace=True)
a2 = result_df.groupby('research_id').sum().reset_index()
a2.to_csv(f"chr{chrom}_count_ancestry_haplotypes_per_person.csv", index=False)
row_sum = a2.iloc[0, 1:].sum()
a2.iloc[:, 1:] = a2.iloc[:, 1:].div(row_sum)
a2.to_csv(f"chr{chrom}_percent_ancestry_haplotypes_per_person.csv", index=False)
