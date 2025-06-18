"""
Plot a heatmap based on phentoype-VNTR association to summarize the entire associations list.

Usage:
python gwas_summary_plot.py --gwas-summary <gwas_input_file>

"""

import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import argparse
from sklearn.metrics import jaccard_score
from scipy.sparse import csr_matrix
from functools import cmp_to_key
from collections import defaultdict
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch



def clustermap(data, col_linkage, filename_base):
    print("Plotting the clustermap")
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

def rename_columns_for_heatmap(data):
    orig_columns = data.columns
    data.columns = data.columns.str.replace('_', ' ')
    data.columns = data.columns.str.lower()
    columns_rename_map = {"skin changes due to solar radiation": "skin changes (solar)",
                          "skin changes due to chronic exposure to nonionizing radiation": "skin changes (NI)",
                          "rheumatic fever and chronic rheumatic heart diseases": "rheumatic heart disease",
                          "hashimoto thyroiditis  chronic lymphocytic thyroiditis ": "hashimoto thyroiditis",
                          }
    data = data.rename(columns=columns_rename_map)
    return data
    

def heatmap(data,
            column_categories,
            filename_base,
            vmax,
            chrom_labels=False,
            category_labels=False,
            xlabel_rotation=45,
            vmin=1e-16):
    print("Plotting the heatmap")
    
    # Locate where each category begins and ends
    xticks, xticklabels = [], []
    axvlines = []
    for idx, category in enumerate(column_categories):
        phenos = column_categories[category]
        min_col = list(data.columns).index(phenos[0])
        max_col = list(data.columns).index(phenos[-1])
        mean_col = int((min_col + max_col) / 2)
        if idx > 0:
            axvlines.append(min_col)
        xticks.append(mean_col)
        xticklabels.append(category)

    # Shorten the name of phenotypes for a better representation in the heatmap
    data = rename_columns_for_heatmap(data)

    plt.clf()
    # Remove annotations due to the large number of rows and columns
    plot = sns.heatmap(data,
                annot=False,
                vmin=vmin,
                vmax=vmax,
                )

    # Separate rows and columns with white lines to distinguish easier
    for idx in range(len(data.columns)):
        plot.axvline(x=idx, color='w',linestyle='-')
    for idx in range(len(data.index)):
        plot.axhline(y=idx, color='w',linestyle='-')

    # Set yticks
    if chrom_labels:
        yticks, yticklabels = [], []
        for chrom in range(1, 23):
            chrom_str = "chr" + str(chrom) + "_"
            # Find the positions for y ticks grouped on chromosomes
            rows_with_chrom = data[data.index.str.startswith(chrom_str)].index
            if len(rows_with_chrom) == 0:
                continue
            min_row = data.index.get_loc(rows_with_chrom[0])
            max_row = data.index.get_loc(rows_with_chrom[-1])
            mean_row = math.ceil((min_row + max_row) / 2)
            yticks.append(mean_row)
            yticklabels.append(chrom_str.replace("_", ""))
        plot.set_yticks(yticks)
        plot.set_yticklabels(yticklabels, fontsize=8)
        # Avoid overlaps by manually moving two closest ones.
        label = plot.get_yticklabels()[-1]
        label.set_position((0, label.get_position()[1] + 10))


    # Place a vertical line at the beginning of each category
    for axvline_pos in axvlines:
        plot.axvline(x=axvline_pos, color='gray',linestyle='--')

    # Set xticks as categories (instead of phenotypes) if indicated
    if category_labels:
        plot.set_xticks(xticks)
        plot.set_xticklabels(xticklabels)
    plt.xticks(rotation=xlabel_rotation, ha="right")

    # Change the automatically generated offset text from e-8 to 10^-8
    cbar = plot.collections[0].colorbar
    text = cbar.ax.yaxis.get_offset_text().get_text()
    text = text.replace("1e", "")
    cbar.ax.yaxis.get_offset_text().set_visible(False)
    cbar.ax.text(-1, 1.04, fr"$\times 10^{{{text}}}$", transform=cbar.ax.transAxes)

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

def get_p_val_df(filename, verbose, default_val, is_continuous=False):
    # Read the summary stats file into a dataframe.
    print("Loading summary statistics into a dataframe.")
    if is_continuous:
        columns = ["Phenotype", "POS", "ID", "REF", "ALT",
                "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ",
                "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"]
    else:
        columns = ["Phenotype", "POS", "ID", "REF", "ALT",
                "PROVISIONAL_REF?", "A1", " OMITTED", "A1_FREQ", "FIRTH?",
                "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT",
                "P", "ERRCODE"]
    summary_df = pd.read_csv(filename, sep="\t", header=0, names=columns)
    summary_df["CHROM"] = summary_df["ID"].apply(lambda x: int(x.split("_")[0].replace("chr", "")))
    summary_df["POS"] = summary_df["POS"].astype(int)
    summary_df["Phenotype"] = summary_df["Phenotype"].apply(lambda x: x.replace("sara_", "")\
                                                                        .replace(".gwas.tab", "")\
                                                                        .replace("_AFR_BLACK", "")\
                                                                        .replace("_AFR_BLACK_female", "")\
                                                                        .replace("_AFR_BLACK_male", "")\
                                                                        .replace("_AMR_HISPANIC", "")\
                                                                        .replace("_AMR_HISPANIC_female", "")\
                                                                        .replace("_AMR_HISPANIC_male", "")\
                                                                        .replace("_EUR_WHITE", "")\
                                                                        .replace("_EUR_WHITE_female", "")\
                                                                        .replace("_EUR_WHITE_male", "")\
                                                                        .replace("_passing_samples_v7.1_gwas.tab", "")\
                                                                        .replace("_passing_samples_v7.1", "")\
                                                                        .split(":")[0]
                                                            )
    phenotypes = summary_df["Phenotype"].unique()
    vntrs = summary_df["ID"].unique()
    vntrs = sorted(list(vntrs), key=cmp_to_key(vntr_order_comp))


    print("Completed loading summary statistics. Now creating p-val df.")
    # Convert it to a dataframe with VNTRs as rows and phenotypes as columns and p-vals as cells.
    # Set the value for insignificant associations to the default value.
    # Defaualt value is 2 * max value, in order to have a mostly bright colored figure.
    print("Current min and max p-values: {}, {}".format(
            summary_df["P"].min(),
            summary_df["P"].max()))
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
    categories = phecode_map_df[['phecode_string', 'phecode_category']]
    phecode_df = phecode_df[phecodes]
    phecode_df = phecode_df.loc[phecodes]
    mapping = dict(zip(phecode_map_df['phecode'], phecode_map_df['phecode_string']))
    max_val = np.max(phecode_df.values)
    phecode_df = phecode_df / max_val
    pheno_df = phecode_df.rename(index=mapping, columns=mapping)
    pheno_df = pheno_df[selected_phenotypes]
    pheno_df = pheno_df.loc[selected_phenotypes]
    for idx, row in pheno_df.iterrows():
        pheno_df.at[idx, idx] = 1
    return pheno_df, categories

def is_sex_specific_pheno(pheno, category):
    pheno = pheno.lower()
    category = category.lower()
    if "_female" in pheno or "_male" in pheno:
        return True
    if pheno.startswith("female") or pheno.startswith("male"):
        return True
    if "pregnancy" in category:
        return True
    if "genitourinary" in category:
        return True
    if "gynecological" in pheno:
        return True
    if "testicular" in pheno or "testis" in pheno or "prostatic" in pheno or "prostate" in pheno:
        return True
    if "uterus" in pheno or "uterine" in pheno or "ovarian" in pheno or "ovary" in pheno:
        return True
    if "vagina" in pheno or "vaginitis" in pheno or "mastopathy" in pheno or "breast" in pheno:
        return True
    if "menstrual" in pheno or "menopausal" in pheno or "menopause" in pheno:
        return True
    if "oligomenorrhea" in pheno or "dysmenorrhea" in pheno or "amenorrhea" in pheno or "endometriosis" in pheno:
        return True
    return False

def distance_based_clustermap_deprecated(pheno_df):
    distance_df = 1 - pheno_df
    condensed_dist = squareform(distance_df.values)
    col_linkage = sch.linkage(condensed_dist, method='average')
    print("col linkage shape: ", np.shape(col_linkage))
    filename_base = "figures/clustermap_{}".format(args.cohort)

    clustermap(data=p_val_df,
        col_linkage=col_linkage,
        filename_base=filename_base)

def filter_vntr_phenos(p_val_df, gwas_phenotypes, categories_map, verbose):
    # Filter VNTRs that appear to have inflated p-values
    # by filtering VNTRs that have very significant p-values in more than a certain number of categories
    p_val_thr_for_filtering = 3E-8
    max_categories_per_vntr = 2
    chosen_vntrs = []
    removed_vntrs = []
    for vntr, row in p_val_df.iterrows():
        categories_w_association = set()
        for pheno in gwas_phenotypes:
            if row[pheno] <= p_val_thr_for_filtering:
                categories_w_association.add(pheno_to_category[pheno])
        if len(categories_w_association) <= max_categories_per_vntr:
            chosen_vntrs.append(vntr)
        else:
            removed_vntrs.append(vntr)
    print("removed_vntrs len {} chosen_vntrs_len {}".format(len(removed_vntrs), len(chosen_vntrs)))
    p_val_df = p_val_df.loc[chosen_vntrs]  
    print("p_val_df shape ", p_val_df.shape)

    # Filter phenotypes that appear to have inflated p-values
    outlier_threshold = 90
    num_vntrs_per_pheno = defaultdict(int)
    for pheno in gwas_phenotypes:
        for vntr in p_val_df.index:
            if p_val_df.loc[vntr, pheno] <= p_val_thr_for_filtering:
                num_vntrs_per_pheno[pheno] += 1
    num_vntrs_threshold = np.percentile(sorted(list(num_vntrs_per_pheno.values())), outlier_threshold)
    if verbose:
        print("num_vntrs_threshold: ", num_vntrs_threshold)
    removed_phenos = {k: v for k, v in num_vntrs_per_pheno.items() if v >= num_vntrs_threshold}
    phenos_to_keep = [item for item in gwas_phenotypes if item not in removed_phenos]
    print("Num phenotypes {} with {} to remove and {} to keep".format(
            len(gwas_phenotypes),
            len(removed_phenos),
            len(phenos_to_keep)))
    p_val_df = p_val_df[phenos_to_keep]

    # Now remove the removed_phenos from the categories map as well
    if verbose:
        print("len categories_map before filtering: ", len(categories_map))
    filtered_categories = defaultdict(list)
    sorted_columns = []
    for category in list(categories_map.keys()):
        for pheno in categories_map[category]:
            if pheno in phenos_to_keep:
                filtered_categories[category].append(pheno)
                sorted_columns.append(pheno)
    categories_map = filtered_categories
    p_val_df = p_val_df[sorted_columns]
    return p_val_df, categories_map

def parse_arguments():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--gwas-summary", required=True,
                         help="Path to the file including the gwas summary results.")
    parser.add_argument("--verbose", action="store_true",
                         help="Enable verbosity.")
    parser.add_argument("--keep-sex-specific-phenos", action="store_true",
                         help="Keep sex specific phenotypes which are otherwise discarded.")
    parser.add_argument("--is-continuous", action="store_true",
                         help="Use this if the phenotype is continuous because the columns " + \
                             "of the significant hits file should be interpretted differently")
    parser.add_argument("--filter", action="store_true",
                         help="Apply filter on VNTRs and phenotypes. " + \
                            "VNTRs that have significant associations in " + \
                            "at least 3 different categories are filtered out. " + \
                            "Also, 10% of phenotypes with most number of significant associations are filtered out.")
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

def update_p_val_categories(p_val_df,
                        phecode_df,
                        gwas_phenotypes,
                        phecode_to_name_file,
                        default_p_val,
                        keep_sex_specific_phenos,
                        verbose):

    pheno_df, category_df = get_pheno_df_from_phecodes(phecode_df=phecode_df,
                                          selected_phenotypes=gwas_phenotypes,
                                          phecode_filename=phecode_to_name_file)

    # Find categories corresponding to each phenotype string
    categories = category_df['phecode_category'].unique()
    sorted_columns = []
    sex_specific_phenotypes = []
    categories_map = defaultdict(list)
    pheno_to_category = {}
    if verbose:
        print("categories head: ", categories[:10])
        print("gwas_phenotypes head: ", gwas_phenotypes[:10])
    for category in categories:
        phecode_strs = category_df[category_df['phecode_category'] == category]["phecode_string"].values
        for phecode_str in phecode_strs:
            # Skip phenotype if we did not find any significant hits
            if phecode_str not in gwas_phenotypes:
                continue
            # Skip sex-specific phenotypes. We often run GWAS on those phenotypes a single sex separately.
            if not keep_sex_specific_phenos and \
                is_sex_specific_pheno(pheno=phecode_str, category=category):
                sex_specific_phenotypes.append(phecode_str)
                continue
            # Skip empty columns
            if p_val_df[phecode_str].min() == default_p_val:
                continue

            categories_map[category].append(phecode_str)
            pheno_to_category[phecode_str] = category
            sorted_columns.append(phecode_str)

    print("p_val_df shape ", p_val_df.shape)
    #print("sex_specific_phenotypes: ", sex_specific_phenotypes)
    print("sorted_columns: ", sorted_columns)
    # Sort columns based on category, for a better appearance in the heatmap and possibly remove some filtered phenotypes.
    p_val_df = p_val_df[sorted_columns]

    # Remove rows (VNTR) with no significant hits after removing the sex-specific phenotypes.
    p_val_df = p_val_df[(p_val_df < default_p_val).any(axis=1)]

    return p_val_df, categories_map

def replace_gene_names(data):
    gene_map = {
        # Significant in AFR
        'chr10_11920339': "UPF2",
        'chr11_33700443': "C11orf91",
        'chr14_44389148': "LINC02307",
        'chr18_47248369': "SKOR2",
        'chr18_48409416': "ZBTB7C",
        'chr20_58981148': "NELFCD",
        'chr21_6103341': "SIK1",
        # Significant in AMR
        'chr1_1449119': "ATAD3C",
        'chr4_3476128': "DOK7",
        'chr6_167668185': "LINC02538",
        'chr7_44036410': "RASA4CP",
        'chr10_59957799': "LINC01553",
        'chr19_40396555': "PRX",
        'chr19_40396820': "PRX",
        'chr19_53538674': "ZNF331",
        # Significant in EUR
        'chr1_165725024': "TMCO1",
        'chr1_165761972': "TMCO1",
        'chr2_21043882': "APOB",
        'chr6_30813385': "LINC00243",
        'chr12_112978257': "OAS2",
        'chr16_89686039': "CDK10",
        'chr16_90007850': "DBNDD1",
        'chr16_90029134': "GAS8",
                }
    return data.rename(index=gene_map)

def main():
    args = parse_arguments()

    # The default value for the matrix for the vntr-phenotype pairs where no significant hit is found.
    #default_val = summary_df["P"].max() * 1.2
    default_p_val = 5e-8
    print("Plotting for cohort ", args.cohort)
    # This is like a placeholder and is larger than the largest observed signficant hit in the matrix.
    p_val_df = get_p_val_df(filename=args.gwas_summary,
                                            default_val=default_p_val,
                                            verbose=args.verbose,
                                            is_continuous=args.is_continuous)
    gwas_phenotypes = list(p_val_df.columns)
    if args.verbose:
        print("gwas_phenotypes head", gwas_phenotypes[:10])
    
    pheno_similarity_file = "pheno_similarities.csv"
    if args.recompute_pheno_similarity:
        phecode_df = get_phenotype_similarity(filename="my_phecode_counts.csv",
                                                  output=pheno_similarity_file)
    else:
        phecode_df = pd.read_csv(pheno_similarity_file, index_col=0)

    # We do not need the categories for blood phenotypes.
    if args.is_continuous:
        categories_map = {}
        # Remove the rows where all values are insignificant.
        p_val_df = p_val_df[(p_val_df < default_p_val).any(axis=1)]
        # A 45 degree rotation is more readable, whereas a 90 degree rotation is ideal for the case of many xticks.
    else:
        # A 45 degree rotation is more readable, whereas a 90 degree rotation is ideal for the case of many xticks.
        p_val_df, categories_map = update_p_val_categories(p_val_df=p_val_df,
                                                        default_p_val=default_p_val,
                                                        phecode_df=phecode_df,
                                                        gwas_phenotypes=gwas_phenotypes,
                                                        phecode_to_name_file=args.phecode_to_name_file,
                                                        keep_sex_specific_phenos=args.keep_sex_specific_phenos,
                                                        verbose=args.verbose)

    
    if args.filter:
        p_val_df, categories_map = filter_vntr_phenos(p_val_df=p_val_df,
                                                      gwas_phenotypes=gwas_phenotypes,
                                                      categories_map=categories_map,
                                                      verbose=args.verbose)
        if args.verbose:
            print("len categories_map after filtering: ", len(categories_map))
            print(categories_map)

    print("p_val_df shape: ", np.shape(p_val_df))
    print("Number of significant values in df: ", (p_val_df < default_p_val).sum().sum())

    filename_base = "figures/heatmap_{}".format(args.cohort)
    
    # Rename VNTR labels (index) with the gene names.
    p_val_df = replace_gene_names(p_val_df)

    heatmap(data=p_val_df,
        column_categories=categories_map,
        filename_base=filename_base,
        vmax=default_p_val)

if __name__ == "__main__":
    main()
