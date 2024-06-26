import pandas as pd
from collections import defaultdict
from cyvcf2 import VCF



def get_hipstr_call():
    # Copying from Nichole's script
    # Copy hipstr calls from the bucket
    # gsutil cp gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/cromwell-execution/targetTR/4519d903-dde6-4c8e-b1ee-bcc2d7cd6dd7/call-sort_index/CBL_test.filtered.sorted.vcf.gz .
    # Extract useful columns from the hipstr call (adjust filenames)
    # !{ echo -en "CHROM\tPOS\tID\t"; bcftools query -l CBL_test.filtered.sorted.vcf.gz | paste -s;bcftools query -f '%CHROM %POS [\t%GB]\n' CBL_test.filtered.sorted.vcf.gz; } > hipstr_CBL.txt
    # Extract useful columns from the imputation call (adjust filenames)
    #bcftools query -r chr11:119206290- -f"%CHROM\t%POS\t%REF\t%ALT\t[%TGT\t]\n" 100_aou_chr11_beagle.vcf.gz|grep -w 119206290 >> genotype_chr11_aou.txt
    call_gt, allele_len
    pass

def read_vcf_file(filename, position):
    print("Working with ", filename)
    variants = []
    reader = VCF(filename)
    samples = reader.samples
    call_gt = {}
    allele_len = defaultdict(list)
    for variant in reader:
        variants.append(variant)
        alt_alleles = variant.ALT
        ref_allele = variant.REF
        if str(variant.POS) != position:
            continue
        print("number of alt alleles ", len(alt_alleles))
        print("len of ref allele", len(ref_allele))
        for i, genotype in enumerate(variant.genotypes):
            if len(genotype) == 2: # No call
                continue
            c1, c2, is_phased = genotype
            for call in [c1, c2]:
                call_gt[samples[i]] = genotype
                if int(call) == 0: #Ref allele
                    allele_len[samples[i]].append(len(ref_allele))
                else:
                    allele_len[samples[i]].append(len(alt_alleles[call - 1]))
                # Sorting the two alleles based on the lengths
                allele_len[samples[i]] = sorted(allele_len[samples[i]])
            #print(genotype, hipstr_allele_len[samples[i]])
    #print("Num variants in the vcf file: ", len(variants))
    return call_gt, allele_len, samples

def get_repeat_count(allele_len, motif_len):
    if motif_len == -1: # In this case, just return allele len
        return allele_len
    # otherwise, return the actual repeat count
    return int(allele_len / motif_len)

def get_diff_motif_count(first_allele_len, second_allele_len, motif_len):
    return abs(
                get_repeat_count(first_allele_len, motif_len) - \
                get_repeat_count(second_allele_len, motif_len) \
                )
            

def compare_pair(first_name, first_allele_len,
                 second_name, second_allele_len,
                 samples, motif_len,
                 repeat_count_threshold,
                 verbose=False):
    is_concordant = []
    for sample in samples:
        # Sanity check that alleles are sorted
        for allele in [first_allele_len[sample], second_allele_len[sample]]:
            if allele[0] > allele[1]:
                print("allele not sorted for sample {}: {} and {}".format(
                    sample,
                    allele[0],
                    allele[1]))
                return

        # Check concordance
        # The difference threshold should hold for both alleles
        if (get_diff_motif_count(first_allele_len[sample][0],
                                second_allele_len[sample][0],
                                motif_len,
                                ) <= repeat_count_threshold) and \
           (get_diff_motif_count(first_allele_len[sample][1],
                                second_allele_len[sample][1],
                                motif_len, 
                                ) <= repeat_count_threshold):
            is_concordant.append(1)
        else:
            if verbose:
                print("example of non-concordant call {} {} and {} {}".format(
                    first_name,
                    first_allele_len[sample],
                    second_name,
                    second_allele_len[sample]))
            is_concordant.append(0)
    concordance = sum(is_concordant)/float(len(samples))

    print("--- Concordance between {} and {} is {:.2f}\%".format(
          first_name, second_name,
          concordance * 100))

    def compare_gt(motif_len, position, repeat_count_thresholds):
    hipstr_filename = "data/CBL_test.filtered.sorted.vcf.gz"
    hipstr_call_gt, hipstr_allele_len, hipstr_samples = read_vcf_file(hipstr_filename, position)
    #caller_filename = "vntr_reference/ref_phased_output_chr15_acan_vntr_apr_22.vcf.gz"
    hipstr_call_gt, hipstr_allele_len, hipstr_samples = read_vcf_file(caller_filename, position)
    
    beagle_filename = "data/CBL_aou_100_phasing_impute_beagle_apr_22_bref.vcf.gz"
    #beagle_filename = "../../../../nichole_imputation/cast-workflows/imputation/wdl/data/output_chr15_10mb_aou_185_lrwgs_w8_o2.vcf.gz"
    #beagle_filename = "../../../../nichole_imputation/cast-workflows/imputation/wdl/data/output_chr15_10mb_aou_275_lrwgs_w8_o2.vcf.gz"
    beagle_call_gt, beagle_allele_len, beagle_samples = read_vcf_file(beagle_filename, position)


    shapeit_filename = "data/CBL_aou_100_phase_shapeit_impute_beagle_apr_22_bref.vcf.gz"
    #shapeit_filename = "data/CBL_aou_100_phase_shapeit_impute_beagle_iflag_apr_22_bref.vcf.gz"
    shapeit_call_gt, shapeit_allele_len, shapeit_samples = read_vcf_file(shapeit_filename, position)

    print("len caller_samples ", len(hipstr_samples))
    print("len beagle_samples ", len(beagle_samples))
    #print("len shapeit_samples ", len(shapeit_samples))
    for repeat_count_threshold in repeat_count_thresholds:
        print("Working with motif length {} and repeat count threshold {}".format(
            motif_len, repeat_count_threshold))
        compare_pair(first_name="caller",
                     first_allele_len=hipstr_allele_len,
                     second_name="beagle",
                     second_allele_len=beagle_allele_len,
                     samples=beagle_samples,
                     motif_len=motif_len,
                     repeat_count_threshold=repeat_count_threshold)
        #compare_pair(first_name="hipstr",
        #             first_allele_len=hipstr_allele_len,
        #             second_name="shapeit",
        #             second_allele_len=shapeit_allele_len,
        #             samples=shapeit_samples)
        #compare_pair(first_name="beagle",
        #             first_allele_len=beagle_allele_len,
        #             second_name="shapeit",
        #             second_allele_len=shapeit_allele_len,
        #             samples=shapeit_samples)

if __name__ == "__main__":
    # For ACAN VNTR
    repeat_count_thresholds = [0, 1, 2, 3]
    position = "88855424"
    motif_len = 57

    # For other test TR
    motif_len = -1
    positions = ["2961502", "8824204", "47310908", "119206290"]
    # For debug
    positions = positions[0:1]
    repeat_count_thresholds = [0]
    for position in positions:
        compare_gt(motif_len, position, repeat_count_thresholds)
