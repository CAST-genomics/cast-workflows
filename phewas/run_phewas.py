from PheTK.Phecode import Phecode
from PheTK.Cohort import Cohort


def call_count_phecodes():
    phecode = Phecode(platform="aou")
    phecode.count_phecode(
        phecode_version="X", 
        icd_version="US",
        phecode_map_file_path=None, 
        output_file_name="my_phecode_counts.csv"
    )


def create_covars(cohort_filename):
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
        output_file_name="cohort_with_covariates.csv"
        )
        
    

def main():
    # The call_count_phecodes only needs to be called once.
    # The phecode_counts csv file is generated for all 250k individuals.
    #call_count_phecodes()
    
    # Before creating covars, a cohort file should be created with the genotype.
    cohort_filename = ""
    create_covars(cohort_filename=cohort_filename)
    

if __name__ == "__main__":
    main()

