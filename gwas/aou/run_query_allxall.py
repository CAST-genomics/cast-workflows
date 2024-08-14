import os
import pandas as pd


def read_region_phenotypes(filename):
    data = pd.read_csv(filename)
    return data

def main():
    data = read_region_phenotypes("region_phenotypes.csv")
    
    # Select a specific region and phenotype pair based on VNTR gene name
    data = data[data["gene"] == "INS"]

    for phenotype in data["phenotypes"][0].split():
        command = "python query_allxall_associations.py " + \
	"--type ACAF " + \
        "--pop META " + \
        "--local " + \
	"--region {} ".format(data["region"][0]) + \
	"--phenoname {} ".format(phenotype) + \
	"--annotate annotation_points.csv"
        os.system(command)

if __name__ == "__main__":
    main()
