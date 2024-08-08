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

def annotate_points(plot, population, region_chrom):
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
                  ]
    # Annotate points
    for point in points:
        chrom, position, label, p_value_dict = point
        # Pick the annotation passed by the user
        if region_chrom != chrom:
            continue
        x = position
        y = p_value_dict[population]
        ax.text(x=x, y=y, s=label, fontsize="medium")
        ax.scatter(x=position, y=y, color="red", marker="^")

def plot_figures(ht, args, region_chrom):
    ## plot a manhattan plot
    filename = "figures/allxall_region_{}_pop_{}_manhattan.png".format(args.phenoname, args.pop)
    title = 'Region P-value for phenotype {}'.format(args.phenoname)
    data = ht.to_pandas()
    plot = sns.relplot(data=data, x="POS", y="Pvalue_log10",
        s=30, aspect=4, linewidth=0, hue="CHR", palette="tab10", legend=None)
    annotate_points(plot, args.pop, region_chrom)
    plot.fig.savefig(filename)

    ## Plot effect sizes
    plt.clf()
    filename = "figures/allxall_{}_pop_{}_effect_sizes.png".format(args.phenoname, args.pop)
    title = 'Region effect sizes for phenotype {}'.format(args.phenoname)
    plot = sns.relplot(data=data, x="POS", y="BETA",
        s=30, aspect=4, linewidth=0, hue="CHR", palette="tab10", legend=None)
    annotate_points(plot, args.pop, region_chrom)
    plot.fig.savefig(filename)


def main():
    args = parse_args()
    #init()

    # Get the corresponding hail path
    ht_path, mt_path = get_hail_path(ex_type=args.type,
                            pop=args.pop,
                            phenoname=args.phenoname)
    ht = hl.read_table(ht_path)
    chrom, start_end = args.region.split(":")
    start, end = start_end.split("-")
    ht = ht.filter(hl.all(
           ht.CHR == chrom,
           ht.POS > int(start),
           ht.POS < int(end)
       ))
    print("len of hail table after filtering: {}".format(ht.count()))
    plot_figures(ht, args, chrom)


if __name__ == "__main__":
    main()

