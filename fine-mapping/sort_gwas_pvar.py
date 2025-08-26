import pandas as pd
import os
import argparse
import sys

"""
Sort GWAS summary statistics to match the order of variants in a PLINK pvar file.

Example command to run:
python sort_gwas_pvar.py --gwas gwas_sumstats.tsv --pvar data.pvar --out output_prefix

"""

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--gwas", help="Name of the gwas summary statistics", type=str, required=True)
    parser.add_argument("--pvar", help="Name of pvar file, ", type=str, required=True)
    parser.add_argument("--out", help="Name of output file", type=str, required=True) 
    args = parser.parse_args()


    # Step 1: Read the pvar file
    pvar = pd.read_csv(args.pvar, sep='\t', comment='#')

    # Step 2: Read the GWAS summary statistics file
    gwas = pd.read_csv(args.gwas, sep='\t')

    gwas['clean_ID'] = gwas['ID']  # default copy

    # Step 3: Handle SNP IDs
    snp_mask = gwas['ID'].str.startswith('SNP:')
    gwas.loc[snp_mask, 'clean_ID'] = (
        gwas.loc[snp_mask, 'ID'].str.replace('SNP:', '', regex=False) + ":" +
        gwas.loc[snp_mask, 'REF'] + ":" +
        gwas.loc[snp_mask, 'ALT']
        )
    gwas['clean_ID'] = gwas['clean_ID'].str.replace(r'\s+', '', regex=True)

    # Step 4: Create a unique identifier for matching
    pvar_id_column = 'ID'  
    gwas_id_column = 'clean_ID' 

    # Step 5: Filter GWAS variants to keep only those in pvar
    common_variants = gwas[gwas[gwas_id_column].isin(pvar[pvar_id_column])]

    # Step 6: Reorder GWAS variants to match pvar order
    common_variants[gwas_id_column] = pd.Categorical(
            common_variants[gwas_id_column],
            categories=pvar[pvar_id_column],
            ordered=True
        )
    
    reordered_gwas= common_variants.sort_values(gwas_id_column)
    reordered_gwas = reordered_gwas.reset_index(drop=True)

    # Step 7: Save the reordered and filtered GWAS summary statistics
    reordered_gwas["clean_ID"].to_csv(f"{args.out}_reordered.txt",index=False,header=None)
    reordered_gwas.to_csv(f"{args.out}_sorted.tsv", sep='\t', index=False)

    print(f"summary for BETA is {reordered_gwas['BETA'].describe()}")
    print(f"summary for SE is {reordered_gwas['SE'].describe()}")
    print(f"Number of variants in pvar: {len(pvar)}")
    print(f"Number of variants in GWAS: {len(gwas)}")
    print(f"Number of common variants: {len(common_variants)}")

if __name__ == "__main__":
    main()