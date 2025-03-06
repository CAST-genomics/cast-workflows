from PheTK.Phecode import Phecode
from PheTK.Cohort import Cohort
from PheTK.PheWAS import PheWAS
from PheTK.Plot import Plot
import gzip
import os
import re
from cyvcf2 import VCF
import pandas as pd
import numpy as np
from collections import Counter
import argparse


def call_count_phecodes():
    phecode = Phecode(platform="aou")
    phecode.count_phecode(
        phecode_version="X", 
        icd_version="US",
        phecode_map_file_path=None, 
        output_file_name="my_phecode_counts.csv"
    )

def _get_alleles(record):
    ref_allele = record.REF
    alt_alleles = record.ALT
    info_fields = record.INFO
    ru = info_fields["RU"]
    if ru is None:
        print("Error: Repeat Unit (RU) not provided in the info field")
        return
    ru_len = len(ru)
    #print("ru_len, ru, len ref allele", ru_len, ru, len(ref_allele))
    len_ref_allele = len(ref_allele)
    len_alt_alleles = [len(alt_allele) for alt_allele in alt_alleles]
    return len_ref_allele, len_alt_alleles, ru_len

def _get_rc_from_allele(allele, ref_allele, alt_alleles, ru_len):
    if allele == 0:
        ref_rc = int(ref_allele/ru_len)
        return ref_rc
    else:
        rc = int(alt_alleles[allele-1]/ru_len)
        return rc

def _write_genotypes_for_locus(tr_vcf, chrom, start, gene, out_filename,
                               vntr_coordinates_threshold=1):
    data = {}
    print("Reading genotypes from the vcf file")
    # Extract the chrom from the tr_vcf file
    words = re.split("_|\.", tr_vcf)
    vcf_chrom = [word for word in words if word.startswith("chr")][0]
    vcf = VCF(tr_vcf)
    vcf_samples = vcf.samples

    all_alleles = []
    empty_calls, no_calls = 0, 0
    variant = list(vcf("{}:{}-{}".format(chrom, int(start)-1, int(start) + 1)))[0]
    print("Processing locus {} {}:{}".format(gene, chrom, start))
    samples_with_calls = set()
    print("Reading ref and alt alleles")
    ref_allele_len, alt_alleles_len, ru_len = _get_alleles(variant)
    print("Reading calls")
    counter = 0
    for i, sample_call in enumerate(variant.genotypes):
        counter += 1
        sample = vcf_samples[i]
        if counter % 10000 == 0:
            print("reading call ", counter)
        rc_1 = _get_rc_from_allele(int(sample_call[0]), ref_allele_len, alt_alleles_len, ru_len)
        rc_2 = _get_rc_from_allele(int(sample_call[1]), ref_allele_len, alt_alleles_len, ru_len)
        rc = (rc_1 + rc_2) / 2.0
        all_alleles.extend([rc_1, rc_2])

        samples_with_calls.add(sample)
        data[sample] = rc
    print("Alleles count for the 4 most common alleles: ", Counter(all_alleles).most_common(4))
    if no_calls + empty_calls > 0:
        print("Skipping {} empty calls and {} no calls for {} on vcf".format(
                empty_calls, no_calls, gene))

    # Write the genotypes in a csv file
    with open(out_filename, "w+") as out_file:
        out_file.write("person_id,genotype\n")
        for person_id in data.keys():
            out_file.write("{},{}\n".format(person_id, data[person_id]))

def create_covars(cohort_filename, output_filename):
    cohort = Cohort(platform="aou", aou_db_version=7)
    cohort.add_covariates(
        cohort_csv_path=cohort_filename,
        natural_age=False,
        age_at_last_event=True,
        sex_at_birth=True,
        ehr_length=False,
        dx_code_occurrence_count=False,
        dx_condition_count=False,
        genetic_ancestry=True,
        first_n_pcs=10,
        drop_nulls=True,
        output_file_name=output_filename
        )
        
    
def run_phewas(locus, cohort_filename, out_filename, min_phecode_count, n_threads, min_cases, verbose):
    covars = ["sex_at_birth",
              "age_at_last_event",
              "pc0", "pc1", "pc2", "pc3", "pc4",
              "pc5", "pc6", "pc7", "pc8", "pc9",
             ]
    genotype_colname = "genotype"
    phewas = PheWAS(
        phecode_version="X",
        phecode_count_csv_path="my_phecode_counts.csv",
        cohort_csv_path=cohort_filename,
        sex_at_birth_col="sex_at_birth",
        male_as_one=True,
        covariate_cols=covars,
        min_phecode_count=min_phecode_count,
        independent_variable_of_interest=genotype_colname,
        min_cases=min_cases,
        output_file_name=out_filename,
        verbose=verbose,
    )
    phewas.run(n_threads=n_threads)

def read_target_loci(filename):
    target_loci = []
    with open(filename, "r") as input_file:
        lines = intput_file.readlines()
    for line in lines:
        chrom, start, gene = line.strip().split(",")
        target_loci.append((str(chrom), str(start), str(gene)))
    return target_loci

def parse_arguments():
    parser = argparse.ArgumentParser(prog = "phewas_runner",
                description="Runs phewas for a given TR locus using pheTK.")

    parser.add_argument("--locus", type=str, required=True,
                        help="The gene name of the locus for phewas.")
    parser.add_argument("--chrom", type=str, required=True,
                        help="The chromosome where the tr is located. " + \
                        "This is used to load the corresponding vcf file to read genotypes.")
    parser.add_argument("--start", type=int, required=True,
                        help="The coordinates where the tr is located." + \
                        "This is used to find the corresponding record on the tr_vcf to load genotypes.")
    parser.add_argument("--n-threads", type=int, default=3,
                        help="Number of threads.")
    parser.add_argument("--min-phecode-count", type=int, default=2,
                        help="Minimum count of a single phecode present for each sample to be considered a case.")
    parser.add_argument("--tr-vcf", type=str, required=True,
                        help="Path to the tr-VCF (imputed) input file.")
    parser.add_argument("--phecode-filename", type=str, default="my_phecode_counts.csv",
                        help="Path to the phecode counts file.")
    parser.add_argument("--outdir", type=str, default="outputs",
                        help="Path to the desired output directory.")
    parser.add_argument("--significance-threshold", type=float, default=5,
                        help="Significance threshold only considered for printing significant hits.")
    parser.add_argument("--min-cases", type=int, default=50,
                        help="Minimum number of cases and controls that should be available " + \
                             "for a phenotype to be included in the phewas.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip the phewas manhattan plot.")
    parser.add_argument("--print-significant-hits", action="store_true",
                        help="Print significant hits.")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Set verbosity to true.")

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    
    # Create output directory if it does not exist
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # The call_count_phecodes only needs to be called once.
    # The phecode_counts csv file is generated for all 250k individuals.
    phecodes_filename = args.phecode_filename
    if not os.path.exists(phecodes_filename):
        call_count_phecodes()

    # Select a locus
    cohort_genotype_filename = "{}/genotypes_{}.csv".format(outdir, args.locus)
    cohort_genotype_covars_filename = "{}/cohort_with_covariates_{}.csv".format(outdir, args.locus)

    # Populate the genotype file
    if os.path.exists(cohort_genotype_filename):
        print("Skipping cohort file building step as the file already exists.")
    else:
        print("Building cohort file.")
        _write_genotypes_for_locus(tr_vcf=args.tr_vcf,
                              chrom=args.chrom,
                              start=args.start,
                              gene=args.locus,
                              out_filename=cohort_genotype_filename)
    
    # Before creating covars, a cohort file should be created with the genotype.
    if os.path.exists(cohort_genotype_covars_filename):
        print("Skipping covars file building step as the file already exists.")
    else:
        print("Building covars file.")
        create_covars(cohort_filename=cohort_genotype_filename,
                      output_filename=cohort_genotype_covars_filename)

    # Run phewas
    phewas_output_filename = "{}/phewas_results_{}_{}.csv".format(
                outdir, args.locus, args.start)
    if os.path.exists(phewas_output_filename):
        print("Skipping run phewas step as the file already exists.")
    else:
        print("Running Phewas.")
        run_phewas(args.locus,
               min_phecode_count=args.min_phecode_count,
               cohort_filename=cohort_genotype_covars_filename,
               out_filename=phewas_output_filename,
               n_threads=args.n_threads,
               verbose=args.verbose,
               min_cases=args.min_cases)

    # Print significant hits.
    if args.print_significant_hits:
        with open(phewas_output_filename, "r") as phewas_result_file:
            for idx, line in enumerate(phewas_result_file.readlines()):
                if idx == 0:
                    # This is the header.
                    continue
                words = line.split(",")
                nlog_pvalue = float(words[4])
                phenotype = words[12]
                category = words[13]
                if float(nlog_pvalue) > args.significance_threshold:
                    print("Significant hit nlog_pvalue {} for VNTR at {}:{} and {} category {}".format(
                        nlog_pvalue,
                        args.locus, args.start,
                        phenotype, category))
                

    # Plot Manhattan
    if not args.no_plot:
        plot_filename = "{}/manhattan_{}_min_phecode_count_{}.png".format(
                    outdir, args.locus, args.min_phecode_count)
        p = Plot(phewas_output_filename)
        p.manhattan(label_values="p_value",
                label_count=10,
                label_value_threshold=2,
                save_plot=True,
                output_file_name=plot_filename,
                output_file_type="png")


if __name__ == "__main__":
    main()



