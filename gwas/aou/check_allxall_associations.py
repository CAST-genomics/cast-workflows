"""
Check if a significant SNP-phenotype association can be found in AoU AllxAll results.
This is meant as a quick check before performing GWAS for TRs, and independent of the TRs.
"""
import os
import argparse
import pandas as pd
#import hail as hl

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
    individual_path = "ht/{ex_type}/{pop}/phenotype_{phenoname}_{ex_type}_results.ht".format(
                        ex_type=ex_type, pop=pop, phenoname=phenoname)
    hail_path = os.path.join(parent_path, individual_path)
    return hail_path

def init():
    hl.init(default_reference = "GRCh38")

def main():
    args = parse_args()
    #init()

    # Get the corresponding hail path
    hail_path = get_hail_path(ex_type=args.type,
                            pop=args.pop,
                            phenoname=args.phenoname)
    print("hail_path: ", hail_path)
    # Read the hail matrix table for the chosen interval
    #mt = hl.read_matrix_table(hail_path)
    ## OR, if complained about matrix table reading a hail table path:
    #ht = hl.read_table(hail_path) and later
    #print(ht.describe())
    #print("len of hail table: {}".format(len(ht)))
    #chrom, start_end = args.region.split(":")
    #start, end = start_end.split("-")
    #ht = ht.filter(hl.all(
    #       ht.CHR=chrom,
    #       ht.POS > int(start),
    #       ht.POS < int(end)
    #   ))
    #print("len of hail table after filtering: {}".format(len(ht)))

    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.region,)])
    #print(mt.describe())
    #print("len of hail matrix table after filtering: {}".format(len(mt)))

    ## plot a manhattan plot
    #filename = "allxall_region_{}_manhattan.pdf".format(args.phenoname)
    #title = 'Region P-value for phenotype {}'.format(args.phenoname)
    #plot = hl.plot.manhattan(mt.info.Pvalue_log10, title=title)
    #plot.fig.savefig(filename)
    ## If save didn't work, try
    #from bokeh.io import save
    # save(plot, filename)
    ## Or try
    #from bokeh.io import export_png
    #export_png(plot, filename=filename)
    # To install bokeh.io, type "pip install bokeh"


if __name__ == "__main__":
    main()

