"""
Check if a significant SNP-phenotype association can be found in AoU AllxAll results.
This is meant as a quick check before performing GWAS for TRs, and independent of the TRs.
"""
import os
import argparse
import pandas as pd
import seaborn as sns
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

def annotate_points(plot):
    ax = plot.ax
    # Annotate points
    for point in [("chr15", 88855424, "ACAN", 41.2),
                  ]:
        chrom, position, label, p_value = point
        x = position
        y = p_value
        ax.text(x=x, y=y, s=label, fontsize="medium")
        ax.scatter(x=position, y=y, color="red", marker="^")

def main():
    args = parse_args()
    #init()

    # Get the corresponding hail path
    ht_path, mt_path = get_hail_path(ex_type=args.type,
                            pop=args.pop,
                            phenoname=args.phenoname)
    print("hail_table_path: ", ht_path)
    print("hail_mt_path: ", mt_path)
    # Read the hail matrix table for the chosen interval
    #mt = hl.read_matrix_table(mt_path)
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.region,)])
    #print(mt.describe())
    #print("len of hail matrix table after filtering: {}".format(len(mt)))
    ## OR, if complained about matrix table reading a hail table path:
    ht = hl.read_table(ht_path)
    print(ht.describe())
    print("len of hail table: {}".format(ht.count()))
    print(ht.show())
    chrom, start_end = args.region.split(":")
    start, end = start_end.split("-")
    ht = ht.filter(hl.all(
           ht.CHR == chrom,
           ht.POS > int(start),
           ht.POS < int(end)
       ))
    print("len of hail table after filtering: {}".format(ht.count()))


    ## plot a manhattan plot
    filename = "allxall_region_{}_manhattan.png".format(args.phenoname)
    title = 'Region P-value for phenotype {}'.format(args.phenoname)
    data = ht.to_pandas()
    print("df ", data)
    print("columns", data.columns)
    plot = sns.relplot(data=data, x="POS", y="Pvalue_log10",
        s=30, aspect=4, linewidth=0, hue="CHR", palette="tab10", legend=None)
    annotate_points(plot)
    plot.fig.savefig(filename)


if __name__ == "__main__":
    main()

