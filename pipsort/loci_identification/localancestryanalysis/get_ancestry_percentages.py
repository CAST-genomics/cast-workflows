import pandas as pd

df = pd.read_csv("chr1_count_ancestry_haplotypes_per_person.csv")
for i in range(2, 23):
    df2 = pd.read_csv(f"chr{i}_count_ancestry_haplotypes_per_person.csv")
    df3 = pd.concat([df, df2])
    df = df3.groupby('research_id').sum().reset_index()

row_sum = df.iloc[0, 1:].sum()
df.iloc[:, 1:] = df.iloc[:, 1:].div(row_sum)
df.to_csv(f"all_ancestry_percentages_per_person.csv", index=False)
