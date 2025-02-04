from PheTK.Phecode import Phecode
from PheTK.Cohort import Cohort
from PheTK.PheWAS import PheWAS
from PheTK.Plot import Plot
import os
import sys
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
import pandas as pd
import argparse

sys.path.append("../utils")
#import utils
from utils import MSG, ERROR

'''
## example code 
!gsutil cp $WORKSPACE_BUCKET/cromwell-execution/targetTR/4519d903-dde6-4c8e-b1ee-bcc2d7cd6dd7/call-sort_index/CBL_test.filtered.sorted.vcf.gz* .
tr_vcf="CBL_test.filtered.sorted.vcf.gz"
python run_pheTK.py --tr-vcf $tr_vcf --region chr11:119206290-119206323 --out CBL_pheTK --n-threads 12

'''
def call_count_phecodes():
    phecode = Phecode(platform="aou")
    phecode.count_phecode(
        phecode_version="X", 
        icd_version="US",
        phecode_map_file_path=None, 
        output_file_name="my_phecode_counts.csv"
    )

#to extract genotype from targetTR (hipstr)
def get_genotype(tr_vcf, region, out_filename):
    # Load TR genotypes
    invcf = utils.LoadSingleReader(tr_vcf, checkgz=True)
    samples = invcf.samples
    region = invcf(region)
    nrecords = 0
    for record in region:
        trrecord = trh.HarmonizeRecord(trh.VcfTypes["hipstr"], record)
        afreqs = trrecord.GetAlleleFreqs()
        genotypes = trrecord.GetLengthGenotypes()
        allele_sum = [sum(item)for item in genotypes]
        trdf = pd.DataFrame({"person_id": samples, "genotype": allele_sum})
        nrecords += 1
    if nrecords == 0:
        ERROR("No matching TR records found")
    if nrecords > 1:
        ERROR("Multiple matching TR records found")

    # Write the genotypes in a csv file
    trdf.to_csv(out_filename, index=False) 
    MSG("Writing out {} genotype file ".format(out_filename))
    

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


def parse_arguments():
    parser = argparse.ArgumentParser(prog = "phewas_runner",
                description="Runs phewas for a given TR locus using pheTK.")

    parser.add_argument("--region",type=str, required=True, help="chr:start-end of target variant")

    parser.add_argument("--n-threads", type=int, default=3,
                        help="Number of threads.")
    parser.add_argument("--min-phecode-count", type=int, default=2,
                        help="Minimum count of a single phecode present for each sample to be considered a case.")
    parser.add_argument("--tr-vcf", type=str, required=True,
                        help="Path to the tr-VCF (imputed) input file.")
    parser.add_argument("--out", type=str, default="stdout", help="Name of output file")
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

   
    cohort_genotype_filename = "{}/genotypes_{}.csv".format(outdir, args.out)
    cohort_genotype_covars_filename = "{}/cohort_with_covariates_{}.csv".format(outdir, args.out)

    # Populate the genotype file
    if os.path.exists(cohort_genotype_filename):
        MSG("Skipping cohort file building step as the file already exists.")
    else:
        MSG("Building cohort file.")
        get_genotype(tr_vcf=args.tr_vcf,
                              region=args.region,
                              out_filename=cohort_genotype_filename)
    
    # Before creating covars, a cohort file should be created with the genotype.
    if os.path.exists(cohort_genotype_covars_filename):
        print("Skipping covars file building step as the file already exists.")
    else:
        print("Building covars file.")
        create_covars(cohort_filename=cohort_genotype_filename,
                      output_filename=cohort_genotype_covars_filename)

    # Run phewas
    phewas_output_filename = "{}/{}_phewas_results_min_phecode_count_{}.csv".format(
                outdir, args.region, args.min_phecode_count)
    if os.path.exists(phewas_output_filename):
        print("Skipping run phewas step as the file already exists.")
    else:
        print("Running Phewas.")
        run_phewas(args.pos,
               min_phecode_count=args.min_phecode_count,
               cohort_filename=cohort_genotype_covars_filename,
               out_filename=phewas_output_filename,
                n_threads=args.n_threads)
    

    # Plot Manhattan
    plot_filename = "{}/manhattan_{}_min_phecode_count_{}.png".format(
                    outdir, args.region, args.min_phecode_count)
    p = Plot(phewas_output_filename)
    p.manhattan(label_values="p_value",
                label_count=10,
                label_value_threshold=2,
                save_plot=True,
                output_file_name=plot_filename,
                output_file_type="png")


if __name__ == "__main__":
    main()



