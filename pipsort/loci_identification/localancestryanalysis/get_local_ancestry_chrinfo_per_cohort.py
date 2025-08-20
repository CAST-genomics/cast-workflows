import numpy as np
import pandas as pd
import sys

chrom = sys.argv[1]
f_cohort = sys.argv[2]
cohortname = f_cohort.split(".")[0]

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


cohort = pd.read_csv(f_cohort)
samples = pd.Series([f"{val}.{i}" for val in cohort['person_id'] for i in range(2)])
df.reset_index(inplace=True)
df.rename(columns={'index':'research_id'}, inplace=True)
samples.name = 'research_id'
samples_df = samples.to_frame()
df = df.merge(samples_df, on='research_id', how='inner')
df.set_index('research_id', inplace=True)

labels={0:'EUR',1:'EAS' ,2:'NAT',3:'AFR', 4:'SAS', 5:'AHG' ,6:'OCE',7:'WAS'}
labels['index']='gnomix_start_pos'

results = {}
for col in df.columns:
    unique_values = np.unique(df[col], return_counts=True)[1]
    results[col] = unique_values

result_df = pd.DataFrame.from_dict(results, orient='index')
result_df.reset_index(inplace=True)
result_df.rename(columns=labels, inplace=True)
result_df.to_csv(f"chr{chrom}_count_ancestry_haplotypes_per_region_{cohortname}.csv", index=False)
row_sum = result_df.iloc[0, 1:].sum()
result_df.iloc[:, 1:] = result_df.iloc[:, 1:].div(row_sum)
result_df.to_csv(f"chr{chrom}_percent_ancestry_haplotypes_per_region_{cohortname}.csv", index=False)

