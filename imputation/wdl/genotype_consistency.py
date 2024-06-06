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
def read_vcf_file(filename):
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

def compare_pair(first_name, first_allele_len, second_name, second_allele_len, samples):
    is_concordant = []
    verbose = False
    for sample in samples:
        if first_allele_len[sample] == second_allele_len[sample]:
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

    print("Concordance between {} and {} is {}".format(first_name, second_name, concordance))

def compare_gt():
    hipstr_filename = "data/CBL_test.filtered.sorted.vcf.gz"
    hipstr_call_gt, hipstr_allele_len, hipstr_samples = read_vcf_file(hipstr_filename)
    
    beagle_filename = "data/CBL_aou_100_phasing_impute_beagle_apr_22_bref.vcf.gz"
    beagle_call_gt, beagle_allele_len, beagle_samples = read_vcf_file(beagle_filename)


    shapeit_filename = "data/CBL_aou_100_phase_shapeit_impute_beagle_apr_22_bref.vcf.gz"
    shapeit_call_gt, shapeit_allele_len, shapeit_samples = read_vcf_file(shapeit_filename)

    print("len hipstr_samples ", len(hipstr_samples))
    print("len beagle_samples ", len(beagle_samples))
    print("len shapeit_samples ", len(shapeit_samples))
    compare_pair(first_name="hipstr",
                 first_allele_len=hipstr_allele_len,
                 second_name="beagle",
                 second_allele_len=beagle_allele_len,
                 samples=beagle_samples)
    compare_pair(first_name="hipstr",
                 first_allele_len=hipstr_allele_len,
                 second_name="shapeit",
                 second_allele_len=shapeit_allele_len,
                 samples=shapeit_samples)
    compare_pair(first_name="beagle",
                 first_allele_len=beagle_allele_len,
                 second_name="shapeit",
                 second_allele_len=shapeit_allele_len,
                 samples=shapeit_samples)

if __name__ == "__main__":
    compare_gt()
