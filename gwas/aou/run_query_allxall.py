import os
import pandas as pd
from time import sleep

def read_region_phenotypes(filename):
    data = pd.read_csv(filename)
    return data

def main():
    data = read_region_phenotypes("region_phenotypes.csv")
   
    chroms = ["chr15:1-102000000", "chr11:1-136000000"]
    for idx, row in data.iterrows():
        print("-------- Running for phenotype ", row["phenotype_name"])
        for chrom in chroms:
            for phenotype in row["phenotypes"].split():
                command = "python query_allxall_associations.py " + \
                "--type ACAF " + \
                "--pop META " + \
                "--region {} ".format(chrom) + \
                "--phenoname {} ".format(phenotype) + \
                "--annotate annotation_points.csv"
                os.system(command)
        sleep(120)

if __name__ == "__main__":
    main()
