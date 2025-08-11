import pandas as pd
import os
import argparse
import sys


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

    # Step 3: Create a unique identifier for matching
    pvar_id_column = 'ID'  
    gwas_id_column = 'clean_ID' 

    # Step 4: Filter GWAS variants to keep only those in pvar
    common_variants = gwas[gwas[gwas_id_column].isin(pvar[pvar_id_column])]

    # Step 5: Reorder GWAS variants to match pvar order
    pvar_order = pvar[pvar_id_column].tolist()
    common_variants[gwas_id_column] = pd.Categorical(
        common_variants[gwas_id_column], categories=pvar_order, ordered=True
    )
    reordered_gwas = common_variants.sort_values(gwas_id_column)

    # Step 6: Save the reordered and filtered GWAS summary statistics
    reordered_gwas["clean_ID"].to_csv(f"{args.out}_reordered.txt",index=False,header=None)
    reordered_gwas.to_csv(f"{args.out}_sorted.tsv", sep='\t', index=False)

    print(f"Number of variants in pvar: {len(pvar)}")
    print(f"Number of variants in GWAS: {len(gwas)}")
    print(f"Number of common variants: {len(common_variants)}")

if __name__ == "__main__":
    main()