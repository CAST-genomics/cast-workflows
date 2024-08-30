import os
import seaborn as sns
import statistics
from matplotlib import pyplot as plt

def filter_vcf(filename, sr_threshold, ml_threshold, out_filename):
    num_updated_calls = 0
    num_calls = 0
    num_prev_no_call = 0
    with open(out_filename, "w+") as output_file:
        with open(filename, "r") as input_file:
            for line in input_file.readlines():
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
                    if word == "./.:0:0:0:0":
                        num_calls += 1
                        num_prev_no_call += 1
                    elif max_likelihood >= ml_threshold and \
                       spanning_read >= sr_threshold:
                        # Passing filters
                        output_words.append(word)
                        num_calls += 1
                    else:
                        # Not passing filters. Replace with a no call.
                        output_words.append("./.:0:0:0:0")
                        num_updated_calls += 1
                        num_calls += 1
                
                # Fix the spacing
                output_line = "\t".join(output_words)
                output_file.write(output_line + '\n')
    print("Filter done with sr_threshold {} ml_threshold {}.".format(
                sr_threshold, ml_threshold))
    print("{} calls were already no call. In addition, {} calls removed ({:.4}% all calls)".format(
                num_prev_no_call,
                num_updated_calls,
                float(num_updated_calls)/num_calls*100,
                ))

if __name__ == "__main__":
    vcf = "data/merged_samples.sorted_batch_1.vcf"
    out_vcf = "data/filtered_merged_samples.sorted_batch_1.vcf"
    for sr_threshold in [4, 6, 8]:
        for ml_threshold in [0.9, 0.95]:
            filter_vcf(filename=vcf,
                   sr_threshold=sr_threshold,
                   ml_threshold=ml_threshold,
                   out_filename=out_vcf)
