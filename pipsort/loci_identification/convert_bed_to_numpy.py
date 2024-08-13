import numpy as np
from pgenlib import PgenReader
import os
import sys
import subprocess



region = sys.argv[1]

bed = region+"_plink.bed"
fam = region+"_plink.fam"
bim = region+"_plink.bim"

num_samples = int(subprocess.check_output(['wc', '-l', fam]).split()[0])
num_variants = int(subprocess.check_output(['wc', '-l', bim]).split()[0])

genotypes = np.empty((num_variants, num_samples), dtype=np.float32)

reader = PgenReader(bytes(str(bed), "utf8"), raw_sample_ct=num_samples, variant_ct=num_variants)

for i in range(0, num_variants):
    reader.read_dosages(i, genotypes[i])

genotypes[genotypes == -9] = 0 #convert missing to reference for now TODO

np.save(region+".npy", genotypes.astype(np.int8))
