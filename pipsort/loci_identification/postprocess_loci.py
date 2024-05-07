import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

chr_regions_file = "chrom_regions_hg38.txt"

with open(chr_regions_file, "r") as f:
    chr_regions = [line.strip() for line in f.readlines()]

chr_regions_dict = {}

for r in chr_regions:
    chrom = r.split(":")[0][3:]
    start, end = r.split(":")[1].split('-')
    chrom = int(chrom)
    start = int(start)
    end = int(end)
    chr_regions_dict[chrom] = (start, end)

loci = pd.read_csv(infile, sep=" ", header=None, names=['chr', 'start', 'end', 'x'])

loci['start'] = loci['start'].apply(lambda x: 1 if x <= 0 else x)

def f(chrom, end):
    actual_end = chr_regions_dict[chrom][1]
    if end > actual_end:
        return actual_end
    else:
        return end

loci['end'] = loci.apply(lambda x: f(x.chr, x.end), axis=1)

loci.to_csv(outfile, sep=" ", header=False, index=False)

