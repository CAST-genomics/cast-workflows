import numpy as np
import pandas as pd
import sys
import subprocess
import os
import matplotlib.pyplot as plt

chrom = sys.argv[1]
bucket = os.getenv('WORKSPACE_BUCKET')

gnomix = f'gnomix-chr{chrom}.msp'
gspath = f"{bucket}/gnomix/outputs/{gnomix}"
if not os.path.exists(gnomix):
    subprocess.run(["gsutil", "cp", gspath, "."], check=True)



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

def prop_different(col):
   arr = col.values.reshape(-1, 2)
   diffs = (arr[:, 0] != arr[:, 1])
   return diffs.mean()

result = df.apply(prop_different, axis=0)
result.to_csv("chr22_prop_het_local_ancestry.txt", header=False)
plt.plot(result.index, result.values)
plt.ylim([0,1])
plt.savefig(f'chr{chrom}_prop_het.png', dpi=300)
