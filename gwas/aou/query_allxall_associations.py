"""
Check if a significant SNP-phenotype association can be found in AoU AllxAll results.
This is meant as a quick check before performing GWAS for TRs, and independent of the TRs.
"""
import os
import argparse
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import hail as hl
import numpy as np

"""
All hail files are located on gs://fc-aou-datasets-controlled/AllxAll/v1/
Hail Matrix Tables: /mt/${POP}_${TYPE}_results.mt
Hail Tables: /ht/${TYPE}/${POP}/phenotype_${PHENONAME}_${TYPE}_results.ht
"""

def parse_args():
    arg_parser = argparse.ArgumentParser(__doc__)

    # Set possible choices for each of the arguments
    ex_type_choices = ["ACAF", "exome", "gene"]
    pop_choices = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "META"]

    # Define arguments
    arg_parser.add_argument("--type", help="Type of the association test",
                            choices=ex_type_choices, type=str, required=True)
    arg_parser.add_argument("--pop", help="Code for the population under study",
                            choices=pop_choices, type=str, required=True)
    arg_parser.add_argument("--phenoname", help="Phenotype under study by concept ID",
                            type=str, required=True)
    arg_parser.add_argument("--region", help="Region to look for SNP-phenotype associations (chr:start-end)",
                            type=str, required=True)
    arg_parser.add_argument("--local", help="Run based on the latest local run, instead of a hail query",
                            action="store_true")
    args = arg_parser.parse_args()
    return args

def get_hail_path(ex_type, pop, phenoname):
    parent_path = "gs://fc-aou-datasets-controlled/AllxAll/v1/"
    individual_ht_path = "ht/{ex_type}/{pop}/phenotype_{phenoname}_{ex_type}_results.ht".format(
                        ex_type=ex_type, pop=pop, phenoname=phenoname)
    hail_table_path = os.path.join(parent_path, individual_ht_path)
    individual_mt_path = "mt/{pop}_{ex_type}_results.mt".format(
                        ex_type=ex_type, pop=pop)
    hail_mt_path = os.path.join(parent_path, individual_mt_path)
    return hail_table_path, hail_mt_path

def init():
    hl.init(default_reference = "GRCh38")

def annotate_points(data, plot, population, region_chrom):
    ax = plot.ax
    points = [ \
            ("chr15", 88855424, "ACAN VNTR", {"AFR": 67.44,
                                         "EUR": 14.45,
                                         "META": 41.2}),
            ("chr11", 2161570, "INS VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr13", 113231980, "CUL4A VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr2", 113130529, "IL1RN VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr12", 56717496, "NACA VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr1", 165761973, "TMCO1 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr8", 116622815, "TMCO1 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr5", 1393581, "SLC6A3-1 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr17", 30221385, "SLC6A4 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr11", 639988, "DRD4 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr5", 132680584, "IL4 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr1", 7829873, "PER3 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr9", 27573484, "C9orf72 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chrX", 71453055, "TAF1 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr21", 43776443 , "CSTB VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr20", 4699397 , "PRNP VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
            ("chr1", 1435799 , "VWA1 VNTR", {"AFR": 0,
                                         "EUR": 0,
                                         "META": 0}),
                  ]
    # Annotate points
    for point in points:
        chrom, position, label, p_value_dict = point
        # Pick the annotation passed by the user
        if region_chrom != chrom:
            continue
        if data["POS"].min() > position:
            continue
        if data["POS"].max() < position:
            continue
        index = len(data[data["POS"] < position])
        x = position
        y = p_value_dict[population]
        ax.text(x=x, y=y, s=label, fontsize="medium")
        ax.scatter(x=x, y=y, color="red", marker="^")

def plot_manhattan(data, x, y, region_chrom,
                   filename, title,
                   args,
                   p_value_threshold=-np.log10(5*10**-8)):

    # Plot a Manhattan plot based on the dataframe
    plot = sns.relplot(data=data, x=x, y=y,
        s=30, aspect=4, linewidth=0,
        hue='CHR', palette="tab10", legend=None,
        )
    chrom_series = data.groupby("CHR")[x].median()

    # Set labels
    plot.ax.set_xlabel("Chromosome")
    plot.ax.set_xticks(chrom_series)
    plot.ax.set_xticklabels(chrom_series.index)
    annotate_points(data, plot, args.pop, region_chrom)

    # Put the threshold
    plot.ax.axhline(p_value_threshold, linestyle="--", linewidth=1)
    plot.fig.savefig(filename, bbox_inches="tight")
    plt.clf()

def plot_figures(data, args, region_chrom):
    #print(data.columns)
    data['ind'] = data.index
    filename = "figures/allxall_manhattan_{}_pop_{}.png".format(args.phenoname, args.pop)
    title = 'Region P-value for phenotype {}'.format(args.phenoname)
    plot_manhattan(data=data,
                   x="POS",
                   y="Pvalue_log10",
                   filename=filename,
                   region_chrom=region_chrom,
                   args=args,
                   title=title)

    ## Plot effect sizes
    filename = "figures/allxall_effect_sizes_{}_pop_{}.png".format(args.phenoname, args.pop)
    title = 'Region effect sizes for phenotype {}'.format(args.phenoname)
    plot_manhattan(data=data,
                   x="POS",
                   y="BETA",
                   filename=filename,
                   region_chrom=region_chrom,
                   args=args,
                   title=title)


def main():
    args = parse_args()
    #init()

    chrom, start_end = args.region.split(":")
    start, end = start_end.split("-")
    start = int(start)
    end = int(end)

    df_dump_filename = "outputs/df_dump_{}_{}.csv".format(chrom, args.phenoname)
    if not args.local:
        # Get the corresponding hail path
        ht_path, mt_path = get_hail_path(ex_type=args.type,
                                pop=args.pop,
                                phenoname=args.phenoname)
        # Read the hail table and filter according to the region
        ht = hl.read_table(ht_path)
        ht = ht.filter(hl.all(
               ht.CHR == chrom,
               ht.POS > int(start),
               ht.POS < int(end)
           ))
        print("len of hail table after filtering: {}".format(ht.count()))
        # Save the dataframe for a later local run.
        data = ht.to_pandas()
        data["POS"] = data["POS"].astype(int)
        data.to_csv(df_dump_filename)
    else:
        # Load data locally
        data = pd.read_csv(df_dump_filename)
        print("# variants loaded", len(data))
        data["POS"] = data["POS"].astype(int)
        data = data[(data["CHR"] == chrom) & \
                    (data["POS"] > start) & \
                    (data["POS"] < end)]
        print("# variants after filtering", len(data))

    plot_figures(data=data, args=args, region_chrom=chrom)


if __name__ == "__main__":
    main()

