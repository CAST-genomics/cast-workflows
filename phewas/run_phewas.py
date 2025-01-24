from PheTK.Phecode import Phecode
from PheTK.Cohort import Cohort
from PheTK.PheWAS import PheWAS
from PheTK.Plot import Plot
import gzip
import os
import pandas as pd
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

def _get_alleles(vcf_df_row):
    ref_allele = vcf_df_row["REF"].array[0]
    alt_alleles = vcf_df_row["ALT"].array[0].split(",")
    info_fields = vcf_df_row["INFO"].array[0].split(";")
    ru = None
    for info_field in info_fields:
        if info_field.startswith("RU"):
            ru = info_field.replace("RU=", "")
    if ru is None:
        print("Error: Repeat Unit (RU) not provided in the info field")
        return
    ru_len = len(ru)
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

def _get_individual_rcs_from_call(call, ref_allele_len, alt_alleles_len, ru_len, sep="/"):
    allele_1, allele_2 = call.split(sep)
    rc_1 = _get_rc_from_allele(int(allele_1), ref_allele_len, alt_alleles_len, ru_len)
    rc_2 = _get_rc_from_allele(int(allele_2), ref_allele_len, alt_alleles_len, ru_len)
    return rc_1, rc_2

def _get_overall_rc_from_call(call, ref_allele_len, alt_alleles_len, ru_len, sep="/"):
    allele_1, allele_2 = call.split(sep)
    rc_1 = _get_rc_from_allele(int(allele_1), ref_allele_len, alt_alleles_len, ru_len)
    rc_2 = _get_rc_from_allele(int(allele_2), ref_allele_len, alt_alleles_len, ru_len)
    return (rc_1 + rc_2)/2.0

def _write_genotypes_for_locus(tr_vcf, chrom, start, gene, out_filename,
                               imputed=True, vntr_coordinates_threshold=10):
    data = {}
    print("Reading genotypes from the vcf file")
    # Read input VCF file into a dataframe
    lines = None

    if tr_vcf.endswith("gz"):
        # It is a compressed vcf file
        with gzip.open(tr_vcf, "rt") as vcf_file:
            lines = vcf_file.readlines()
    elif tr_vcf.endswith("vcf"):
        # It is an uncompressed vcf file
        with open(tr_vcf, "r") as vcf_file:
            lines = vcf_file.readlines()
    else:
        print("Error: Cannot recognize tr-vcf file format. Should be either a vcf file or a vcf.gz file")
        exit(1)

    # Remove header lines for the dataframe
    lines_clean = [line for line in lines if not line.startswith("##")]
    columns = lines_clean[0].replace("\n", "").replace("#", "").split("\t")
    df_lines = [line.split("\t") for line in lines_clean[1:]]
    vcf_df = pd.DataFrame(df_lines, columns=columns, dtype=str)
    vcf_df["CHROM"] = vcf_df["CHROM"].astype(str)
    vcf_df["POS"] = vcf_df["POS"].astype(str)
    print("DF shape", vcf_df.shape)

    # Get all calls for this locus
    all_alleles = []
    empty_calls, no_calls = 0, 0
    locus_calls = vcf_df[(vcf_df["CHROM"] == chrom) & \
                         (vcf_df["POS"].astype(int) < start + vntr_coordinates_threshold) & \
                         (vcf_df["POS"].astype(int) > start - vntr_coordinates_threshold)]
    if len(locus_calls) == 0:
        # In a different chromosome
        print("Error: Locus in gene {} not found in {}.".format(gene, chrom))
        return
    elif len(locus_calls) > 1:
        print("Error: Multiple loci within the {}bp of input coordinates. ".format(
                    vntr_coordinates_threshold) + \
                "Please adjust the vntr_coordinates_threshold to have one " + \
                "locus selected.")
        return
    print("Processing annotation {} {}:{}".format(gene, chrom, start))
    samples_with_calls = set()
    shared_columns = []
    print("Reading ref and alt alleles")
    ref_allele_len, alt_alleles_len, ru_len = _get_alleles(locus_calls)
    print("Reading calls")
    counter = 0

    for column in vcf_df.columns:
        counter += 1
        if counter % 10000 == 0:
            print("reading call ", counter)
        if column.isnumeric():
            if imputed:
                if len(locus_calls[column]) == 0:
                    # Equal to no call
                    continue
                call = locus_calls[column].array[0]
                call = call.split(":")[0]
                #print("Locus call after split ", call)
                rc = _get_overall_rc_from_call(call, ref_allele_len, alt_alleles_len, ru_len, sep="|")
                # Get individual alleles for plotting
                rc_1, rc_2 = _get_individual_rcs_from_call(call, ref_allele_len, alt_alleles_len, ru_len, sep="|")
                all_alleles.extend([rc_1, rc_2])
            else: # not imputed
                if len(locus_calls[column].array) == 0:
                    # An error causing an empty value in the vcf file
                    empty_calls += 1
                    continue
                call = locus_calls[column].array[0]
                call = call.split(":")[0]
                if call == ".":
                    # Equal to no call
                    no_calls += 1
                    continue
                rc = _get_overall_rc_from_call(call, ref_allele_len, alt_alleles_len, ru_len, sep="/")
                # Get individual alleles for plotting
                rc_1, rc_2 = _get_individual_rcs_from_call(call, ref_allele_len, alt_alleles_len, ru_len, sep="/")
                all_alleles.extend([rc_1, rc_2])


            samples_with_calls.add(column)
            data[column] = rc
        else:
            shared_columns.append(column)
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
        
    
def run_phewas(locus, cohort_filename, out_filename, min_phecode_count, n_threads):
    covars = ["sex_at_birth",
              #"genetic_ancestry",
              "age_at_last_event",
              "pc0", "pc1", "pc2", "pc3", "pc4",
              "pc5", "pc6", "pc7", "pc8", "pc9",
             ]
    genotype_colname = "genotype"
    example_phewas = PheWAS(
        phecode_version="X",
        phecode_count_csv_path="my_phecode_counts.csv",
        cohort_csv_path=cohort_filename,
        sex_at_birth_col="sex_at_birth",
        male_as_one=True,
        covariate_cols=covars,
        min_phecode_count=min_phecode_count,
        independent_variable_of_interest=genotype_colname,
        min_cases=50,
        output_file_name=out_filename,
    )
    example_phewas.run(n_threads=n_threads)

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
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    # The call_count_phecodes only needs to be called once.
    # The phecode_counts csv file is generated for all 250k individuals.
    #call_count_phecodes()
    
    # Create output directory if it does not exist
    outdir = "outputs"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Select a locus
    cohort_genotype_filename = "{}/genotypes_{}.csv".format(outdir, args.locus)
    cohort_genotype_covars_filename = "{}/cohort_with_covariates_{}.csv".format(outdir, args.locus)

    # Populate the genotype file
    tr_vcf = "../imputation/wdl/data/imputed_{}.annotated.rh.vcf.gz".format(args.chrom)
    if not os.path.exists(cohort_genotype_filename):
        _write_genotypes_for_locus(tr_vcf=tr_vcf,
                              chrom=args.chrom,
                              start=args.start,
                              gene=args.locus,
                              out_filename=cohort_genotype_filename,
                              imputed=True)
    
    # Before creating covars, a cohort file should be created with the genotype.
    if not os.path.exists(cohort_genotype_covars_filename):
        create_covars(cohort_filename=cohort_genotype_filename,
                      output_filename=cohort_genotype_covars_filename)

    # Run phewas
    phewas_output_filename = "{}/{}_phewas_results_min_phecode_count_{}.csv".format(
                outdir, args.locus, args.min_phecode_count)
    if not os.path.exists(phewas_output_filename):
        run_phewas(args.locus,
               min_phecode_count=args.min_phecode_count,
               cohort_filename=cohort_genotype_covars_filename,
               out_filename=phewas_output_filename,
                n_threads=args.n_threads)
    

    # Plot Manhattan
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



