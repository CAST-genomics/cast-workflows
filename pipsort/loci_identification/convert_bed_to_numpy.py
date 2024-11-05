import numpy as np
from pgenlib import PgenReader
import os
import sys
import subprocess



bed = sys.argv[1]
num_samples = int(sys.argv[2])
num_variants = int(sys.argv[3])

fileprefix = bed[:-4] #chop off .bed

genotypes = np.empty((num_variants, num_samples), dtype=np.float32)

reader = PgenReader(bytes(str(bed), "utf8"), raw_sample_ct=num_samples, variant_ct=num_variants)

for i in range(0, num_variants):
    reader.read_dosages(i, genotypes[i])

#genotypes[genotypes == -9] = 0 #convert missing to reference for now TODO

np.save(fileprefix+".npy", genotypes.astype(np.int8))
