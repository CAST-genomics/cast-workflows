"""
Plot a heatmap based on phentoype-VNTR association to summarize the entire associations list.

Usage:
python gwas_summary_plot.py --gwas-summary <gwas_input_file>

"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
from functools import cmp_to_key
from collections import defaultdict



def heatmap(data, filename_base):
    print("Plotting the heatmap")
    plt.clf()
    # Remove annotations due to the large number of rows and columns
    plot = sns.heatmap(data,
                annot=False,
                xticklabels=False,
                yticklabels=False)
    yticks, yticklabels = [], []
    # Find the positions for y ticks grouped on chromosomes
    for chrom in range(1, 23):
        chrom_str = "chr" + str(chrom) + "_"
        rows_with_chrom = data[data.index.str.startswith(chrom_str)].index
        if len(rows_with_chrom) == 0:
            continue
        min_row = data.index.get_loc(rows_with_chrom[0])
        max_row = data.index.get_loc(rows_with_chrom[-1])
        mean_row = int((min_row + max_row) / 2)
        yticks.append(mean_row)
        yticklabels.append(chrom_str.replace("_", ""))
        #print("row number for the first value ", data.index.get_loc(rows_with_chrom[0]))
        #print("row number for the last value ", data.index.get_loc(rows_with_chrom[-1]))
        #print(f"min and max row for chrom {chrom}: {min(rows_with_chrom)}, {max(rows_with_chrom)}")
    plot.set_yticks(yticks, yticklabels)

    plt.savefig(filename_base + ".pdf", bbox_inches="tight")
    plt.savefig(filename_base + ".png", bbox_inches="tight")

def vntr_order_comp(a, b):
    # Returns -1 if a < b, 1 if a > b and 0 if a = b
    chr_a, pos_a = a.split("_")
    chr_b, pos_b = b.split("_")
    chr_a = int(chr_a.replace("chr", ""))
    chr_b = int(chr_b.replace("chr", ""))
    pos_a = int(pos_a)
    pos_b = int(pos_b)
    # First compare the chromosomes
    if chr_a < chr_b:
        return -1
    if chr_a > chr_b:
        return 1
    # Then compare the positions within the same chromosome
    if pos_a < pos_b:
        return -1
    if pos_a > pos_b:
        return 1
    return 0

def get_p_val_df(filename, verbose):
    # Read the summary stats file into a dataframe.
    print("Loading summary statistics into a dataframe.")
    columns = ["Phenotype", "POS", "ID", "REF", "ALT",
                "PROVISIONAL_REF?", "A1", " OMITTED", "A1_FREQ", "FIRTH?",
                "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT",
                "P", "ERRCODE"]
    summary_df = pd.read_csv(filename, sep="\t", header=0, names=columns)
    summary_df["CHROM"] = summary_df["ID"].apply(lambda x: int(x.split("_")[0].replace("chr", "")))
    summary_df["POS"] = summary_df["POS"].astype(int)
    summary_df["Phenotype"] = summary_df["Phenotype"].apply(lambda x: x.replace("sara_", "")\
                                                                        .replace("_EUR_WHITE.gwas.tab", "")\
                                                                        .replace("_AFR_BLACK.gwas.tab", "")\
                                                                        .replace("_passing_samples_v7.1_gwas.tab", "")\
                                                                        .split(":")[0]
                                                            )
    phenotypes = summary_df["Phenotype"].unique()
    vntrs = summary_df["ID"].unique()
    vntrs = sorted(list(vntrs), key=cmp_to_key(vntr_order_comp))


    print("Completed loading summary statistics. Now creating p-val df.")
    # Convert it to a dataframe with VNTRs as rows and phenotypes as columns and p-vals as cells.
    # Set the value for insignificant associations to the default value.
    # Defaualt value is 2 * max value, in order to have a mostly bright colored figure.
    default_val = summary_df["P"].max() * 1.2
    p_val_df = pd.DataFrame(columns=phenotypes, index=vntrs, dtype=float).fillna(default_val)


    multiple_entries = defaultdict(list)
    index = 0
    for vntr in vntrs:
        index += 1
        for phenotype in phenotypes:
            p_val = summary_df[(summary_df["ID"] == vntr) & \
                                                        (summary_df["Phenotype"] == phenotype)]["P"]
            # No significant association between the VNTR and phenotype.
            # Skin and keep the default val.
            if isinstance(p_val, pd.core.series.Series) and len(p_val) == 0:
                continue
            # More than one entry for the VNTR-phenotype pair
            elif isinstance(p_val, pd.core.series.Series) and len(p_val) > 1:
                multiple_entries[phenotype].append((vntr, tuple(list(p_val))))
                p_val = p_val.head(1)
            p_val_df.loc[vntr, phenotype] = p_val.item()
    print("P-val dataframe populated.")

    if verbose:
        print(f"More than one entry for {len(multiple_entries)} phenotypes. Picking the first one.")
        print(f"full list {multiple_entries}" )

    return p_val_df

def get_phenotype_similarity(filename):
    print("Reading the phenotype counts file")
    pheno_df = pd.read_csv(filename)
    print(pheno_df)

    samples = pheno_df["person_id"].unique()
    phecodes = pheno_df["phecode"].unique()
    print("samples len ", len(samples))
    print("phecode len ", len(phecodes))

    print("Creating a case-control dataframe")
    case_control_df = pd.DataFrame(index=samples, columns=phecodes)
    index = 0
    for sample in samples:
        if index >= 10:
            break
        for phecode in phecodes:
            index += 1
            count = pheno_df[(pheno_df["person_id"] == sample) & \
                             (pheno_df["phecode"] == phecode)]["count"]
            if isinstance(count, pd.core.series.Series) and len(count) == 0:
                # No entry for this sample and this phecode. So it's a control
                case_control_df.loc[sample, phecode] = 0
            else:
                if count >= 2:
                    case_control_df.loc[sample, phecode] = 1
                else:
                    case_control_df.loc[sample, phecode] = 0
    print(case_control_df)            
                
            

def parse_arguments():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--gwas-summary", required=True,
                         help="Path to the file including the gwas summary results.")
    parser.add_argument("--verbose", action="store_true",
                         help="Enable verbosity.")
    parser.add_argument("--cohort", required=True,
                         help="Cohort for which the gwas summary stats are provided. " +
                                "Choose among ALL, EUR, AFR, AMR")

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    
    phenotypes_df = get_phenotype_similarity("my_phecode_counts.csv")
    p_val_df = get_p_val_df(filename=args.gwas_summary,
                            verbose=args.verbose)


    filename_base = "heatmap_{}".format(args.cohort)
    heatmap(data=p_val_df,
            filename_base=filename_base)

if __name__ == "__main__":
    main()
