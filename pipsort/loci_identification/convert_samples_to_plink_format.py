import pandas as pd
import sys

samples_file = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(samples_file)
samples = pd.DataFrame(0, index=df.index, columns=['FID'])
samples['IID'] = df['person_id']
samples.to_csv(outfile, sep="\t", index=False, header=False)
