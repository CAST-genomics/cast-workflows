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
import warnings

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
    arg_parser.add_argument("--type", help="Type of the association test.",
                            choices=ex_type_choices, type=str, required=True)
    arg_parser.add_argument("--pop", help="Code for the population under study.",
                            choices=pop_choices, type=str, required=True)
    arg_parser.add_argument("--phenoname", help="Phenotype under study by concept ID.",
                            type=str, required=True)
    arg_parser.add_argument("--region", help="Region to look for SNP-phenotype associations (chr:start-end). " + \
                            "Does not currently provide whole genome.",
                            type=str, required=True)
    arg_parser.add_argument("--annotate", help="csv file with points to use as annotations on the plots.",
                            type=str, required=False)
    arg_parser.add_argument("--local", help="Run based on the latest local run, instead of a hail query.",
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

def annotate_points(data, plot, population, region_chrom, annotate_filename):
    ax = plot.ax
    points_df = pd.read_csv(annotate_filename)
    # Annotate points
    for idx, point in points_df.iterrows():
        position = point["position"]
        # Pick the annotation passed by the user within the chosen region
        if region_chrom is not None:
            if region_chrom != point["chrom"]:
                continue
            if data["POS"].min() > position:
                continue
            if data["POS"].max() < position:
                continue
        index = len(data[data["POS"] < position])
        x = position
        y = point[population]
        ax.text(x=x, y=y, s=point["name"], fontsize="medium")
        ax.scatter(x=x, y=y, color="red", marker="^")

def plot_manhattan(data, x, y,
                   region_chrom,
                   filename, title,
                   pop,
                   annotate_filename,
                   p_value_threshold=-np.log10(5*10**-8)):
    # Filter warning for tight layout in matplotlib
    warnings.filterwarnings("ignore", module="seaborn\..*")

    # Plot a Manhattan plot based on the dataframe
    plot = sns.relplot(data=data, x=x, y=y,
        s=30, aspect=4, linewidth=0,
        hue='CHR', palette="tab10", legend=None,
        )

    # Set labels
    if region_chrom == None:
        # Plot is for the whole chromosome, so we need xticks to
        # refer to chromosomes rather than coordinates
        plot.ax.set_xlabel("Chromosome")
        chrom_series = data.groupby("CHR")[x].median()
        plot.ax.set_xticks(chrom_series)
        plot.ax.set_xticklabels(chrom_series.index)
    else:
        # Plot is for a region in one specific chromosome.
        plot.ax.set_xlabel(region_chrom)

    # Annotate VNTR positions
    if annotate_filename is not None:
        annotate_points(data, plot, pop, region_chrom, annotate_filename)

    # Put the threshold
    plot.ax.axhline(p_value_threshold, linestyle="--", linewidth=1)
    plot.fig.savefig(filename, bbox_inches="tight")
    plt.clf()

def plot_figures(data, pop, phenoname, annotate_filename, region_chrom):
    #print(data.columns)
    data['ind'] = data.index
    filename = "figures/allxall_manhattan_{}_pop_{}.png".format(phenoname, pop)
    title = 'Region P-value for phenotype {}'.format(phenoname)
    plot_manhattan(data=data,
                   x="POS",
                   y="Pvalue_log10",
                   filename=filename,
                   region_chrom=region_chrom,
                   pop=pop,
                   annotate_filename=annotate_filename,
                   title=title,)
    output_filename = "outputs/allxall_dataframe_{}_{}_pop_{}.csv".format(region_chrom, phenoname, pop)
    data.to_csv(output_filename)

def main():
    args = parse_args()

    if args.region is not None:
        # Plotting for the selected region
        chrom, start_end = args.region.split(":")
        start, end = start_end.split("-")
        start = int(start)
        end = int(end)
        df_dump_filename = "outputs/df_dump_{}_{}.csv".format(chrom, args.phenoname)
    else:
        # Plotting for the whole genome
        chrom, start, end = None, None, None
        df_dump_filename = "outputs/df_dump_wg_{}.csv".format(args.phenoname)

    if not os.path.exists("outputs"):
        os.mkdir("outputs")
    if not os.path.exists("figures"):
        os.mkdir("figures")

    
    if not args.local:
        init()
        # Get the corresponding hail path
        ht_path, mt_path = get_hail_path(ex_type=args.type,
                                pop=args.pop,
                                phenoname=args.phenoname)
        # Read the hail table and filter according to the region
        ht = hl.read_table(ht_path)
        if chrom is not None:
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
        print("Num variants loaded: ", len(data))
        data["POS"] = data["POS"].astype(int)
        if chrom is not None:
            data = data[(data["CHR"] == chrom) & \
                    (data["POS"] > start) & \
                    (data["POS"] < end)]
        print("Num variants after filtering for position:", len(data))

    plot_figures(data=data,
                 pop=args.pop,
                 phenoname=args.phenoname,
                 annotate_filename=args.annotate,
                 region_chrom=chrom)


if __name__ == "__main__":
    main()

