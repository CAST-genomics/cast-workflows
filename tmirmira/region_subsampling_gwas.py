import numpy as np
import os
import pandas as pd
from pgenlib import PgenReader
from sklearn.preprocessing import normalize
import subprocess
import sys



plink_prefix = sys.argv[1]


if not os.path.exists(f"{plink_prefix}.npy"):
    bed = f"{plink_prefix}.bed"
    bim = f"{plink_prefix}.bim"
    fam = f"{plink_prefix}.fam"
    result = subprocess.run(["wc", "-l", bim], capture_output=True, text=True, check=True)
    num_variants = int(result.stdout.strip().split()[0])
    result = subprocess.run(["wc", "-l", fam], capture_output=True, text=True, check=True)
    num_samples = int(result.stdout.strip().split()[0])

    genotypes = np.empty((num_variants, num_samples), dtype=np.float32)

    reader = PgenReader(bytes(str(bed), "utf8"), raw_sample_ct=num_samples, variant_ct=num_variants)

    for i in range(0, num_variants):
        reader.read_dosages(i, genotypes[i])

    #genotypes[genotypes == -9] = np.nan #encode missing
    genotypes[genotypes == -9] = 0 #convert missing to reference for now TODO
    np.save(f"{plink_prefix}.npy", genotypes.astype(np.int8))


start = time.time()
m = np.load(f"{plink_prefix}.npy")
m = m.T #need to transpose because it will be numvariants x numsamples and we want the transpoe

print(m.shape)

m = normalize(m, axis=0).astype(np.float16)
print(m.dtype)
end = time.time()
print("time to load and preprocess matrix = ", end - start)
