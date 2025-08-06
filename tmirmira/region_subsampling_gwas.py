import gc
import numpy as np
import os
import pandas as pd
from pgenlib import PgenReader
from sklearn.preprocessing import normalize
import subprocess
import sys
import time
from sklearn.preprocessing import StandardScaler



plink_prefix = sys.argv[1]


bed = f"{plink_prefix}.bed"
bim = f"{plink_prefix}.bim"
fam = f"{plink_prefix}.fam"
result = subprocess.run(["wc", "-l", bim], capture_output=True, text=True, check=True)
num_variants = int(result.stdout.strip().split()[0])
print("num variants: ", num_variants)
result = subprocess.run(["wc", "-l", fam], capture_output=True, text=True, check=True)
num_samples = int(result.stdout.strip().split()[0])
print("num samples: ", num_samples)

if not os.path.exists(f"{plink_prefix}.bin"):

    #genotypes = np.empty((num_variants, num_samples), dtype=np.float32)


    genotypes = np.memmap(f"{plink_prefix}.bin", dtype=np.int8, mode='w+', shape=(num_variants, num_samples))

    reader = PgenReader(bytes(str(bed), "utf8"), raw_sample_ct=num_samples, variant_ct=num_variants)

    tmp = np.empty(num_samples, dtype=np.float32)
    for i in range(num_variants):
        reader.read_dosages(i, tmp)
        tmp[tmp == -9] = 0
        genotypes[i] = tmp.astype(np.int8)

    genotypes.flush()
    del genotypes
    gc.collect()
    #genotypes[genotypes == -9] = np.nan #encode missing
    #genotypes[genotypes == -9] = 0 #convert missing to reference for now TODO
    #np.save(f"{plink_prefix}.npy", genotypes.astype(np.int8))

    print("done writing file")


start = time.time()
genotypes = np.memmap(f"{plink_prefix}.bin", dtype=np.int8, mode='r', shape=(num_variants, num_samples))
# If you want a full array in RAM:
m = np.array(genotypes, dtype=np.float32)
m = m.T #need to transpose because it will be numvariants x numsamples and we want the transpoe

print(m.shape)

scaler = StandardScaler(with_mean=False)
m = scaler.fit_transform(m).astype(np.float16)

#m = normalize(m, axis=0).astype(np.float16)
print(m.dtype)
end = time.time()
print("time to load and preprocess matrix = ", end - start)
