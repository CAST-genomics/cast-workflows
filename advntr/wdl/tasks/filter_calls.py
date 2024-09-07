import os
import seaborn as sns
import pandas as pd
import statistics
from collections import defaultdict
from matplotlib import pyplot as plt

def filter_vcf(filename, sr_threshold, ml_threshold, out_filename, verbose=False):
    num_updated_calls = 0
    num_calls = 0
    num_prev_no_call = 0
    num_unique_alleles = []
    with open(out_filename, "w+") as output_file:
        with open(filename, "r") as input_file:
            for line in input_file.readlines():
                unique_alleles = defaultdict(int)
                # Copy header lines as is
                if line.startswith("#"):
                    output_file.write(line)
                    continue
                
                # Variant line
                output_words = []
                words = line.strip().split()
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
                    else:
                        # Not passing filters. Replace with a no call.
                        output_words.append("./.:0:0:0:0")
                        num_updated_calls += 1
                        num_calls += 1
                
                # Fix the spacing
                output_line = "\t".join(output_words)
                output_file.write(output_line + '\n')
                # Update unique alleles
                num_unique_alleles.append(len(unique_alleles))
    if verbose:
        print("Filter done with sr_threshold {} ml_threshold {}.".format(
                sr_threshold, ml_threshold))
        print("Total {} calls. {} calls were already no call. In addition, {} calls removed ({:.4}% all calls)".format(
                num_calls,
                num_prev_no_call,
                num_updated_calls,
                float(num_updated_calls)/num_calls*100,
                ))
    return float(num_updated_calls)/num_calls*100, num_unique_alleles

def plot_heatmap(data, filename):
    plt.clf()
    ax = sns.heatmap(data=data, annot=True)
    ax.set_title("Percent calls removed")
    ax.get_figure().savefig(filename)

def filter_search(vcf, out_vcf):
    results_df = pd.DataFrame(columns=["sr_threshold", "ml_threshold", "percent removed"])
    for sr_threshold in range(12):
        for ml_threshold in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1]:
            percent_removed, _ = filter_vcf(filename=vcf,
                   sr_threshold=sr_threshold,
                   ml_threshold=ml_threshold,
                   out_filename=out_vcf)
            results_df.loc[len(results_df)] = [sr_threshold, ml_threshold, percent_removed]
    results_df = results_df.pivot(
         index="sr_threshold",
         columns="ml_threshold",
         values="percent removed",
            )
    #print(results_df)
    plot_heatmap(data=results_df,
         filename="filter_heatmap.png")

def plot_hist(data, threshold, filename):
    plt.clf()
    bins = list(range(threshold + 1))
    bins_shifted = [pos-0.5 for pos in bins[1:]]
    bins_labels = [str(pos) for pos in bins[1:]]
    #ax = sns.histplot(data,
    #                  stat="percent",
    #plt.hist(data,
    #                  bins=bins,
    ax = sns.ecdfplot(data,
                      stat="proportion",
                      complementary=True)
    #ax = plt.gca()
    ax.set_xticks(bins_shifted, bins_labels)
    ax.set_title("Allele distribution capped at {}".format(threshold))
    ax.set(xlabel="Number of unique alleles", ylabel="Proportion of VNTRs")
    ax.get_figure().savefig(filename)

if __name__ == "__main__":
    #vcf = "data/merged_samples.sorted_batch_1.vcf"
    vcf = "data/merged_samples_p_g_vntrs_all.sorted.vcf"
    #out_vcf = "data/filtered_merged_samples.sorted_batch_1.vcf"
    out_vcf = "data/filtered_merged_samples_p_g_vntrs_all.sorted.vcf"
    #filter_search(vcf, out_vcf)
    _, alleles_dist = filter_vcf(filename=vcf,
                   sr_threshold=0,
                   ml_threshold=0,
                   out_filename=out_vcf)
    threshold = 20
    alleles_dist_capped = [min(threshold, item) for item in alleles_dist]
    plot_hist(alleles_dist_capped, threshold, "allele_distribution.png")

