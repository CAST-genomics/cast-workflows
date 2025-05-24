"""
Plot a heatmap based on phentoype-VNTR association to summarize the entire associations list.

Usage:
python gwas_summary_plot.py --gwas-summary <gwas_input_file>

"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
from sklearn.metrics import jaccard_score
from scipy.sparse import csr_matrix
from functools import cmp_to_key
from collections import defaultdict
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch





def heatmap(data, col_linkage, filename_base):
    print("Plotting the heatmap")
    plt.clf()
    # Remove annotations due to the large number of rows and columns
    plot = sns.clustermap(data,
                col_linkage=col_linkage,
                annot=False,
                row_cluster=False,
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
    plot.ax_heatmap.set_yticks(yticks)
    plot.ax_heatmap.set_yticklabels(yticklabels)
    #plot.ax.ax_heatmap.axes.yticks(yticks, yticklabels)

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

def get_phenotype_similarity(filename, output):
    print("Reading the phenotype counts file")
    pheno_df = pd.read_csv(filename)

    samples = pheno_df["person_id"].unique()
    phecodes = pheno_df["phecode"].unique()
    print("samples len ", len(samples))
    print("phecode len ", len(phecodes))

    # Create a case-control dataframe.
    print("Creating a case-control dataframe")
    case_control_df = pd.DataFrame(index=samples, columns=phecodes, dtype=bool)
    case_control_df[:] = False
    num_cases = defaultdict(int)
    index = 0
    for idx, row in pheno_df.iterrows():
        phecode = row["phecode"]
        #if phecode not in ["272.1", "EM_239.1", "3027114", "MS_718", "EM_239.11"]:
        #    continue
        #print(f"processing phecode {phecode}")
        if idx % 100000 == 0:
            print(f"processing row {idx}")
        sample = row["person_id"]
        count = row["count"]
        if count >= 2:
            case_control_df.loc[sample, phecode] = True
            num_cases[phecode] += 1


    # Find the distance between two phenotypes
    # Convert boolean DataFrame to a sparse matrix for memory and compute efficiency
    X = csr_matrix(case_control_df.values.astype(np.uint8))
    
    # Compute the intersection: dot product gives number of co-True values per pair
    intersection = np.dot(X.T, X)  # shape: (4000, 4000)
    print("intersection shape: ", intersection.shape)
    
    # Convert diagonal to a dense array (number of True values per column)
    col_sums = intersection.diagonal()
    print("col_sums shape: ", col_sums.shape)
    
    # Compute union: |A ∪ B| = |A| + |B| - |A ∩ B|
    union = col_sums[:, None] + col_sums[None, :] - intersection.toarray()
    print("union shape: ", union.shape)
    
    # Compute Jaccard similarity: |A ∩ B| / |A ∪ B|
    jaccard_similarity = np.where(union != 0, intersection.toarray() / union, 0.0)
    #jaccard_similarity = intersection.toarray() / union
    
    print("jaccard_similarity type", type(jaccard_similarity))
    print("jaccard_similarity shape", jaccard_similarity.shape)
    similarity_df = pd.DataFrame(jaccard_similarity, columns=phecodes, index=phecodes)
    similarity_df.to_csv("pheno_similarities.csv")
    return similarity_df
            
def get_pheno_df_from_phecodes(phecode_df, selected_phenotypes,
                                          phecode_filename):
    phecode_map_df = pd.read_csv(phecode_filename, header=0)
    phecodes = phecode_map_df['phecode']
    print(phecode_map_df)
    phecode_df = phecode_df[phecodes]
    phecode_df = phecode_df.loc[phecodes]
    mapping = dict(zip(phecode_map_df['phecode'], phecode_map_df['phecode_string']))
    max_val = np.max(phecode_df.values)
    phecode_df = phecode_df / max_val
    pheno_df = phecode_df.rename(index=mapping, columns=mapping)
    pheno_df = pheno_df[selected_phenotypes]
    pheno_df = pheno_df.loc[selected_phenotypes]
    print(pheno_df)
    for idx, row in pheno_df.iterrows():
        pheno_df.at[idx, idx] = 1
    return pheno_df

def parse_arguments():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--gwas-summary", required=True,
                         help="Path to the file including the gwas summary results.")
    parser.add_argument("--verbose", action="store_true",
                         help="Enable verbosity.")
    parser.add_argument("--recompute-pheno-similarity", action="store_true",
                         help="Recompute phenotype similarities based on a Jaccard similarity of case-controls. " + \
                            "Warning: This could take several hours.")
    parser.add_argument("--cohort", required=True,
                         help="Cohort for which the gwas summary stats are provided. " +
                                "Choose among ALL, EUR, AFR, AMR")
    parser.add_argument("--phecode-to-name-file", required=True,
                         help="Path to the file containing phecodes and phenotype names")

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    p_val_df = get_p_val_df(filename=args.gwas_summary,
                            verbose=args.verbose)
    gwas_phenotypes = list(p_val_df.columns)
    print("gwas_phenotypes", gwas_phenotypes)
    
    pheno_similarity_file = "pheno_similarities.csv"
    if args.recompute_pheno_similarity:
        phecode_df = get_phenotype_similarity(filename="my_phecode_counts.csv",
                                                  output=pheno_similarity_file)
    else:
        phecode_df = pd.read_csv(pheno_similarity_file, index_col=0)
        print(phecode_df)
    pheno_df = get_pheno_df_from_phecodes(phecode_df=phecode_df,
                                          selected_phenotypes=gwas_phenotypes,
                                          phecode_filename=args.phecode_to_name_file)
    distance_df = 1 - pheno_df
    condensed_dist = squareform(distance_df.values)
    col_linkage = sch.linkage(condensed_dist, method='average')
    print("col linkage shape: ", np.shape(col_linkage))

    
    print("p_val_df shape: ", np.shape(p_val_df))

    filename_base = "heatmap_{}".format(args.cohort)
    heatmap(data=p_val_df,
            col_linkage=col_linkage,
            filename_base=filename_base)

if __name__ == "__main__":
    main()
