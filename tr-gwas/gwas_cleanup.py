#!/usr/bin/env python3
"""
Script to clean up GWAS summary statistics
 

./gwas_cleanup.py \
--phenotype hdl_cholesterol_AMR_HISPANIC_gwas.tab

"""

import argparse
import json
import os
import subprocess
import sys
import pandas as pd
import numpy as np
from google.cloud import storage

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
		path = os.getenv("WORKSPACE_BUCKET")+"/tr-gwas_result/"+item.strip()
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
    df['P'] = df['P'].replace('NA', np.nan)
    df_cleaned = df.dropna(subset=['P'])
    print(f'{df_cleaned.shape[0]} number of STRs are left after removing NA')
    phenoname = file_path.replace(".tab", ".tsv")
    df_cleaned.to_csv(f"{outdir}/{phenoname}", sep='\t', index=False)
    
    return phenoname

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
    outdir = "outputs"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Read data
    if args.phenotype is not None:
        phenotype_list = GetGWASPath(args.phenotype)
    else:
        phenotype_list  = [blob.name for blob in bucket_name.list_blobs(prefix="tr-gwas_result/") if blob.name.endswith('.tab')]
    for phenotype in phenotype_list:
        pheno_name = phenotype.split("/")[-1]
        if os.path.exists(pheno_name):
             sys.stdout.write(f"{pheno_name} already exists locally")
        else:
            DownloadGWAS(phenotype)
           
        tsv_file = Cleanupfile(pheno_name,outdir)
    # Compress and index   
        CompressIndex(tsv_file,outdir)

if __name__ == "__main__":
    main()
