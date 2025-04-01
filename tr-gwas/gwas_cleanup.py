#!/usr/bin/env python3
"""
Script to clean up GWAS summary statistics
 

./gwas_cleanup.py 

"""

import argparse
import json
import os
import subprocess
import sys
import pandas as pd
from google.cloud import storage

sys.path.append("../utils")
import aou_utils

def GetGWASPath(phenotype):
	"""
	Download a GCP path locally
	Arguments
	---------
	phenotype : str
	   GCP path

	"""
	phenotype_array = []
	phenotypes = [item.strip() for item in phenotype.split(',')]
	for item in phenotypes: 
		path = os.getenv("WORKSPACE_BUCKET")+"/tr-gwas_result/"+item.strip()+"_gwas.tab"
		phenotype_array.append(path)	
	return phenotype_array

def DownloadGWAS(file_path): 
    """
    Download a GCP path locally

    Arguments
    ---------
    filename : str
       GCP path
    """
    cmd = f"gsutil cp {file_path} ."
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    print(output.decode("utf-8"))   
    

def Cleanupfile(file_path,outdir):  
    """
	Read GWAS summary statistics file and remove extra header
	Arguments
	---------
	phenotype : str
	   GCP path

	"""
    data = []
    # Initialize a flag for reading data after the header
    is_reading_data = False
    columns = []
    header_found = False  # To track if we've already found the first header
    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    # Open the file
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespace
            if line.startswith('==>'):
                continue
            if line.startswith('#CHROM'):
                if not header_found:
                    columns = line.split()  # Store the first header
                    header_found = True  # Mark that we've found the first header
                continue  # Skip all subsequent header lines

            if header_found:
                # Only process non-empty lines
                if line:
                    data.append(line.split())  # Split the line by whitespace or tab
    df = pd.DataFrame(data, columns=columns)
    df = df[df['TEST'] == 'ADD']
    df.to_csv(f"{outdir}/{file_path}_tsv", sep='\t', index=False)
    return (f"{file_path}_tsv")

def CompressIndex(file,outdir):  
    cmd = f"bgzip {outdir}/{file}"
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    print(output.decode("utf-8"))
    cmd = f"tabix -p vcf {outdir}/{file}.gz"
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    print(output.decode("utf-8"))
              

def main():
    parser = argparse.ArgumentParser(description="Remove extra header and compress and index gwas summary statistics.")
    parser.add_argument("--phenotype", help="name of the phenotype, seperated by comma", required=False, type=str)
    args = parser.parse_args()

    bucket_name = os.getenv("WORKSPACE_BUCKET")
    # Create output directory if it does not exist
    outdir = "gwas/outputs"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Read data
    if args.phenotype is not None:
        phenotype_list = GetGWASPath(args.phenotype)
    else:
        phenotype_list  = [blob.name for blob in bucket_name.list_blobs(prefix="tr-gwas_result/") if blob.name.endswith('.tab')]
    for phenotype in phenotype_list:
        if not os.path.exists(phenotype):
            DownloadGWAS(phenotype)
            phenotype = phenotype.split("/")[-1]
        
        tsv_file = Cleanupfile(phenotype,outdir)
    # Compress and index   
        CompressIndex(tsv_file,outdir)

if __name__ == "__main__":
    main()
