import sys
from pgenlib import PgenWriter
from pgenlib import PgenReader
import numpy as np
import pandas as pd


snps_infile = sys.argv[1]
strs_infile = sys.argv[2]
snp_pvar_file = sys.argv[3]
str_pvar_file = sys.argv[4]
snp_ct = int(sys.argv[5])
str_ct = int(sys.argv[6])
variant_ct = snp_ct + str_ct
outpgen = sys.argv[7]
outpvar = sys.argv[8]
num_samples = int(sys.argv[9])


pgen_writer = PgenWriter(bytes(str(outpgen), "utf8"), 
        num_samples, 
        variant_ct=variant_ct, 
        dosage_present=True)


dosages = np.empty(num_samples, dtype=np.float32)



snp_pvar = pd.read_csv(snp_pvar_file, sep="\t")
str_pvar = pd.read_csv(str_pvar_file, sep="\t")
snp_reader = PgenReader(bytes(str(snps_infile), "utf8"))
str_reader = PgenReader(bytes(str(strs_infile), "utf8"))
snp_i = 0
str_i = 0

outpvar_f = open(outpvar, 'w')
header_str = '\t'.join(snp_pvar.columns)+"\n"
outpvar_f.write(header_str)


for i in range(0, variant_ct):
    curr_snp_loc = np.inf
    curr_str_loc = np.inf
    if snp_i < snp_ct:
        curr_snp_loc = snp_pvar.iloc[snp_i]['POS']
    if str_i < str_ct:
        curr_str_loc = str_pvar.iloc[str_i]['POS']
    if curr_snp_loc <= curr_str_loc:
        snp_reader.read_dosages(snp_i, dosages)
        outpvar_f.write("\t".join(snp_pvar.iloc[snp_i,:].apply(str).values)+"\n")
        snp_i += 1
    else:
        str_reader.read_dosages(str_i, dosages)
        outpvar_f.write("\t".join(str_pvar.iloc[str_i,:].apply(str).values)+"\n")
        str_i += 1
    pgen_writer.append_dosages(dosages)

outpvar_f.close()
snp_reader.close()
str_reader.close()
pgen_writer.close()