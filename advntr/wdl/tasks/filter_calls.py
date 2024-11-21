import os
import gzip
import sys
import seaborn as sns
import pandas as pd
import numpy as np
import statistics
from scipy.stats import ks_2samp
from collections import defaultdict
from matplotlib import pyplot as plt

def read_vcf(filename, locus):
    alleles_distribution_dict = defaultdict(int)#lambda: sys.float_info.epsilon)
    alleles_distribution = []
    num_no_call = 0
    if filename.endswith("gz"):
        input_file = gzip.open(filename, "rt")
        lines = input_file.readlines()
    else:
        input_file = open(filename, "r")
        lines = input_file.readlines()

    for line_idx, line in enumerate(lines):
        if line.startswith("#"):
            continue
        
        # Variant line
        output_words = []
        words = line.strip().split()
        has_at_least_one_call = False
        if locus is not None:
            locus_chr, locus_pos = locus
            if words[0].strip() != locus_chr or \
               words[1].strip() != str(locus_pos):
                   continue
        for idx, word in enumerate(words):
            if idx < 9:
                continue
            genotype = word.split(":")[0]
            if "/" in genotype:
                sep = "/"
            elif "|" in genotype:
                sep = "|"
            else:
                print("Genotype {} does not have valid seperator".format(genotype))
            alleles = word.split(":")[0].split(sep)
            if alleles[0] == "." or alleles[1] == ".":
                #print("call is ", word)
                num_no_call += 1
                pass
            else:
                alleles_distribution_dict[int(alleles[0])] += 1
                alleles_distribution_dict[int(alleles[1])] += 1
    max_allele_rc = int(max(alleles_distribution_dict.keys()))
    print("num_no_call is {}".format(num_no_call))
    for idx in range(max_allele_rc):
        alleles_distribution.append(float(alleles_distribution_dict[idx]))
    return alleles_distribution

def filter_vcf(filename, sr_threshold, ml_threshold, out_filename, verbose=False):
    num_updated_calls = 0
    num_calls = 0
    num_prev_no_call = 0
    num_all_no_calls = 0
    num_vntrs_in_chrs = defaultdict(int)
    num_unique_alleles = []
    polymorphic_vntrs = []
    if filename.endswith("gz"):
        input_file = gzip.open(filename, "rt")
        lines = input_file.readlines()
    else:
        input_file = open(filename, "r")
        lines = input_file.readlines()
    with open(out_filename, "w+") as output_file:
        for line in lines:
            unique_alleles = defaultdict(int)
            # Copy header lines as is
            if line.startswith("#"):
                output_file.write(line)
                continue
            
            # Variant line
            output_words = []
            words = line.strip().split()
            has_at_least_one_call = False
            for idx, word in enumerate(words):
                if idx < 9:
                    output_words.append(word)
                    continue
                max_likelihood = float(word.split(":")[4])
                spanning_read = int(word.split(":")[2])
                if word.startswith("./."):
                    num_calls += 1
                    num_prev_no_call += 1
                    num_updated_calls += 1
                    output_words.append(word)
                elif max_likelihood >= ml_threshold and \
                   spanning_read >= sr_threshold:
                    # Passing filters
                    output_words.append(word)
                    num_calls += 1
                    alleles = word.split(":")[0].split("/")
                    if alleles[0] == "." or alleles[1] == ".":
                        print("call is ", word)
                    else:
                        unique_alleles[int(alleles[0])] += 1
                        unique_alleles[int(alleles[1])] += 1
                    has_at_least_one_call = True
                else:
                    # Not passing filters. Replace with a no call.
                    output_words.append("./.:0:0:0:0")
                    num_updated_calls += 1
                    num_calls += 1
            if not has_at_least_one_call:
                num_all_no_calls += 1

            # Fix the spacing
            output_line = "\t".join(output_words)
            output_file.write(output_line + '\n')

            # Update unique alleles
            num_unique_alleles.append(len(unique_alleles))
            if len(unique_alleles) > 1:
                polymorphic_vntrs.append(words[2])
            if words[0] == "chr15" and words[1] == "88855424":
                print("For ACAN, num unique alleles: ", len(unique_alleles))
            if len(unique_alleles) >= 100 and verbose:
                print("For {}:{} we have {} unique alleles with min {} max {}".format(
                        words[0], words[1], len(unique_alleles),
                        min(unique_alleles),
                        #np.mean(unique_alleles),
                        max(unique_alleles)))
            num_vntrs_in_chrs[words[0]] += 1

    # Write polymorphic VNTR IDs
    with open("polymorphic_vntrs.txt", "w+") as vntr_id_file:
        for vntr in polymorphic_vntrs:
            vntr_id_file.write(vntr + "\n")
    if verbose:
        print("Filter done with sr_threshold {} ml_threshold {}.".format(
                sr_threshold, ml_threshold))
        print("Total {} calls. {} calls were already no call. In addition, {} calls removed ({:.4}% all calls)".format(
                num_calls,
                num_prev_no_call,
                num_updated_calls,
                float(num_updated_calls)/num_calls*100,
                ))
        print("{} loci had all none / no calls".format(num_all_no_calls))
        for key in num_vntrs_in_chrs.keys():
            print("Total number of VNTRs in {}: {}".format(key, num_vntrs_in_chrs[key]))
    return float(num_updated_calls)/num_calls*100, num_unique_alleles

def plot_heatmap(data, filename):
    plt.clf()
    ax = sns.heatmap(data=data, annot=True)
    ax.set_title("Percent calls removed")
    ax.get_figure().savefig(filename)

def filter_search(vcf, out_vcf):
    results_df = pd.DataFrame(columns=["sr_threshold", "ml_threshold", "percent removed"])
    #for sr_threshold in range(12):
    #    for ml_threshold in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1]:
    for sr_threshold in [6]:
        for ml_threshold in [0.95]:
            percent_removed, alleles_dist = filter_vcf(filename=vcf,
                   sr_threshold=sr_threshold,
                   ml_threshold=ml_threshold,
                   out_filename=out_vcf,
                    verbose=True)
            results_df.loc[len(results_df)] = [sr_threshold, ml_threshold, percent_removed]
    results_df = results_df.pivot(
         index="sr_threshold",
         columns="ml_threshold",
         values="percent removed",
            )
    #print(results_df)
    threshold=100
    num_h_calls = sum([1 for item in alleles_dist if item == 1])
    print("num loci with only 1 allele: ", num_h_calls)
    num_no_calls = sum([1 for item in alleles_dist if item == 0])
    print("num loci with no call: ", num_no_calls)
    num_too_many_calls = sum([1 for item in alleles_dist if item >= threshold])
    print("num loci with >= {} alleles: {}".format(threshold, num_too_many_calls))
    print([item for item in alleles_dist if item >= threshold])
    print("max unique loci:" , max(alleles_dist))
    plot_heatmap(data=results_df,
         filename="filter_heatmap.png")

def plot_hist(data, threshold, filename, is_cumulative, is_reverse):
    plt.clf()
    bins = list(range(threshold + 1))
    #print(bins)
    bins_shifted = [pos + 0.5 for pos in bins]
    bins_labels = [str(pos) for pos in bins]
    if not is_cumulative:
        ax = sns.histplot(data,
                      bins=bins,
                      stat="density",
                      )
        ax.set_xticks(bins_shifted, bins_labels)
    elif not is_reverse:
        ax = sns.histplot(data,
                      bins=bins,
                      stat="density",
                      cumulative=True,
                      )
        ax.set_xticks(bins_shifted, bins_labels)
    else:
        ax = sns.ecdfplot(data,
                      stat="proportion",
                      complementary=True)
    ax.set_title("Allele distribution capped at {}".format(threshold))
    ax.set(xlabel="Number of unique alleles", ylabel="Proportion of VNTRs")
    ax.get_figure().savefig(filename)

def plot_allele_histogram(vcf, out_vcf):
    _, alleles_dist = filter_vcf(filename=vcf,
                   sr_threshold=6,
                   ml_threshold=0.95,
                   out_filename=out_vcf,
                   verbose=True)
    threshold = 20
    num_h_calls = sum([1 for item in alleles_dist if item == 1])
    print("num loci with only 1 allele: ", num_h_calls)
    num_no_calls = sum([1 for item in alleles_dist if item == 0])
    print("num loci with no call: ", num_no_calls)
    num_too_many_calls = sum([1 for item in alleles_dist if item >= threshold])
    print("num loci with >= {} alleles: {}".format(threshold, num_too_many_calls))
    #print("max uniq alleles: ", max(alleles_dist))
    alleles_dist_capped = [min(threshold, item) for item in alleles_dist]
    plot_hist(alleles_dist_capped, threshold, "unique_allele_distribution.png",
                is_cumulative=False, is_reverse=False)
    plot_hist(alleles_dist_capped, threshold, "unique_allele_distribution_cumulative.png",
                is_cumulative=True, is_reverse=False)
    plot_hist(alleles_dist_capped, threshold, "unique_allele_distribution_reverse_cumulative.png",
                is_cumulative=True, is_reverse=True)

def validate_genotypes(vcf1, vcf2, verbose=False):
    alleles_dist1 = read_vcf(filename=vcf1,
                   locus=("chr15", "88855424"),
                   )
    alleles_dist2 = read_vcf(filename=vcf2,
                   locus=("chr15", "88855424"),
                   )
    if False:
        # extend the distributions with 0 values
        if len(alleles_dist1) < len(alleles_dist2):
            len_diff = len(alleles_dist2) - len(alleles_dist1)
            alleles_dist1.extend([sys.float_info.epsilon] * len_diff)
        elif len(alleles_dist2) < len(alleles_dist1):
            len_diff = len(alleles_dist1) - len(alleles_dist2)
            alleles_dist2.extend([sys.float_info.epsilon] * len_diff)
        if len(alleles_dist1) != len(alleles_dist2):
            print("Error! two allele distributions are not the same length after padding.")
    if verbose:
        print("for dist 1 len {} min {} mean {} max {}".format(
            len(alleles_dist1),
            min(alleles_dist1),
            statistics.mean(alleles_dist1),
            max(alleles_dist1)))
        print("for dist 2 len {} min {} mean {} max {}".format(
            len(alleles_dist2),
            min(alleles_dist2),
            statistics.mean(alleles_dist2),
            max(alleles_dist2)))
    diff = ks_2samp(alleles_dist1, alleles_dist2, method="exact")
    print("ks_2samp between the two distributions is: {}".format(diff))


if __name__ == "__main__":
    vcf = "data/merged_all_lrwgs_p_g_vntrs.sorted.vcf.gz"
    vcf_prev_ref = "data/ACAN_merged_samples.sorted.vcf"
    vcf_srwgs = "data/merged_samples_all_rh.sorted.vcf.gz"
    out_vcf = "data/filtered_sr_6_ml_95_merged_all_lrwgs_p_g_vntrs_v2.sorted.vcf"
    #out_vcf = "data/tmp_filtered_merged_all_p_g_vntrs.sorted.vcf"
    #vcf="data/filtered_merged_all_p_g_vntrs.sorted.vcf"
    filter_search(vcf, out_vcf)
    #plot_allele_histogram(vcf, out_vcf)
    #validate_genotypes(vcf1=vcf, vcf2=vcf_prev_ref, verbose=True)
