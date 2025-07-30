import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import subprocess

bucket = os.getenv('WORKSPACE_BUCKET')

chrom = sys.argv[1]

centromeres = pd.read_csv("centromere_positions.csv", index_col='chrom')
centromere_start = int(centromeres.loc[f'chr{chrom}']['chromStart'])
centromere_end = int(centromeres.loc[f'chr{chrom}']['chromEnd'])

gnomix = f'gnomix-chr{chrom}.msp'

if not os.path.exists(gnomix):
    subprocess.run(
        f"gsutil cp {bucket}/gnomix/outputs/{gnomix} ./",
        shell=True,
        check=True
    )

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

df.columns = df.columns.astype(int)  # or float if needed
ancestry_counts = df.apply(pd.Series.value_counts, axis=0).fillna(0)
ancestry_counts = ancestry_counts.loc[[0, 1, 2, 3, 4, 7]]
ancestry_percentages = ancestry_counts / df.shape[0]

# Plot
plt.figure(figsize=(18, 6))

labels = {0: 'EUR', 1: 'EAS', 2: 'NAT', 3: 'AFR', 4: 'SAS', 7: 'WAS'}
colors = {0: 'dodgerblue', 1: 'orange', 2: 'green', 3: 'red', 4: 'purple', 7: 'brown'}

for anc in ancestry_percentages.index:
    plt.plot(
        ancestry_percentages.columns,
        ancestry_percentages.loc[anc],
        label=labels[anc],
        color=colors[anc]
    )

plt.xlabel("Position")
plt.ylabel("Proportion")
plt.title(f"Chromosome {chrom}")
plt.axvline(x=centromere_start, color='red', linestyle='--', linewidth=1)
plt.axvline(x=centromere_end, color='red', linestyle='--', linewidth=1)
plt.legend()
plt.tight_layout()
plt.savefig(f"overall_ancestry_proportion_chrom{chrom}.png")

