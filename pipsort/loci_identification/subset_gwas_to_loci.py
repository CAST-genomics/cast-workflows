import sys
import pandas as pd
import os

chrom = sys.argv[1]
from_bp = float(sys.argv[2])
to_bp = float(sys.argv[3])
gwas_file = sys.argv[4]
subset_gwas_file = sys.argv[5]

gwas = pd.read_csv(gwas_file, sep="\t", comment='#')
gwas = gwas[(gwas['chrom']=="chr"+chrom) & (gwas['pos']>=from_bp) & (gwas['pos']<=to_bp)]
gwas.to_csv(subset_gwas_file, sep="\t", index=False)
