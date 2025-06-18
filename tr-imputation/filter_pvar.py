import re
from collections import defaultdict
import argparse
import sys


def parse_dsfield(ds_string):
    try:
        alleles_part, counts_part = ds_string.split(";")
        alleles = [float(x) for x in alleles_part.split(",")]
        counts= [int(x) for x in counts_part.split(",")]
        if len(alleles) != len(counts):
            raise ValueError("DSCOUNT format error")

        d = defaultdict(float)
        for allele,count in zip(alleles,counts):
            d[allele] += count
            total_count = sum(counts)
        if total_count == 0:
            raise ValueError("Total count is zero")
            
        for allele in d:
            d[allele] /= total_count
        return d
        
    except Exception as e:
        print(f"DSCOUNT parse error: {ds_string} → {e}")
        return {}
    

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--input", help="Input file name", required=True, type=str)
    parser.add_argument("--output", help="Output file name", required=True, type=str)
    parser.add_argument("--MAF", help="Minimun of minor allels frequency", type=float, default=0.01)
    parser.add_argument("--min-alleles", help="Minimun of unique allele numbers", type=int, default=2)

    args = parser.parse_args()

    passing_ids = []

    with open(args.input, "r") as fin:
        for line in fin:
            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.strip().split()
            id_field = fields[2]
            info_field = fields[7:]
            # Extract DSCOUNT field
            # info_fields is a list of one item, so take that string
            info_string = info_field[0]  
            #print(info_string)
            info_parts = info_string.split("=")[2:][0]
    
            ds_counts = parse_dsfield(info_parts)
            #print(ds_counts)
            passing_allele_count = sum(1 for count in ds_counts.values() if count >= args.MAF)
            if passing_allele_count >= args.min_alleles:
                passing_ids.append(id_field)


    # Write passing IDs to output
    with open(args.output, "w") as fout:
        for id in passing_ids:
            fout.write(f"{id}\n")

    print(f" Wrote {len(passing_ids)} passing locus IDs to: {args.output}")

if __name__ == "__main__":
	main()
