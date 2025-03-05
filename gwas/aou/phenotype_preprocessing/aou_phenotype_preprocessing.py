#!/usr/bin/env python3

"""
Preprocess AoU phenotypes.
Output cleaned phenotype and covariates files
for the specified phenotype

Requires these environment variables to be set:
WORKSPACE_CDR
WORKSPACE_BUCKET
"""

import argparse
import aou_queries
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import datetime

# To avoid the SettingWithCopyWarning
pd.options.mode.chained_assignment = None

SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")

def MSG(msg_str):
    """
    Write a helpful progress message to stderr

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    sys.stderr.write("[aou_phenotype_preprocessing]: %s\n"%msg_str.strip())

def ERROR(msg_str):
    """
    Write an error message to stderr then quit
    with exit code 1

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    MSG("ERROR: " + msg_str)
    sys.exit(1)

def SQLToDF(sql):
    """
    Extract Google bigquery results to pandas df

    Arguments
    ---------
    sql : str
       Query string

    Returns
    -------
    df : pandas.DataFrame
       Dataframe with the results
    """
    df = pd.read_gbq(
        sql,
        dialect="standard",
        use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
        progress_bar_type="tqdm_notebook")
    return df

def OverlapDrugMeasurement(measurement_datetime, start, end):
    """
    Determine whether a phenotype measurement overlaps
    the start/end time of a drug

    Currently: counts as no if phenotype measurement
    is before drug start

    Arguments
    ---------
    measurement_datetime : datetime
       Time phenotype measurement taken
    start : datetime
       Time drug started
    end : datetime
       Time drug ended

    Returns
    -------
    res : int
       0 (not on drug) or 1 (on drug)
    """
    if measurement_datetime < start:
        res = 0
    else:
        res = 1
    return res

def GetDrugData(concept_id):
    """
    Extract drug data for a concept id

    Arguments
    ---------
    concept_id : int
       AoU concept ID

    Returns
    -------
    drugdata : pandas.DataFrame
       Includes columns: person_id, start (drug start), end (drug end)
    """
    drug_sql = aou_queries.ConstructDrugExposureSQL(concept_id)
    drugdata = SQLToDF(drug_sql)
    drugdata = drugdata.groupby(["person_id"]).agg(start=('drug_exposure_start_datetime', np.min), \
        end=('drug_exposure_end_datetime', np.max)).reset_index()
    return drugdata


def my_median(series):
    """
    compute median of median 
    odd number use np.median
    even number pick either one

    Arguments
    ---------
    series : pandas.DataFrame column


    Returns
    -------
    median: pandas.DataFrame
       
    """
    if len(series)%2==1:
        return np.median(series)
    else:
        my_median = sorted(series)[int(len(series)/2)]
        return my_median


def parse_phecode_cohort_file(samples, demog, args):
    # Create a samples and demog df
    all_samples_data = pd.merge(samples, demog, on="person_id", how="left")

    # Set the age to the current age
    all_samples_data["age"] = datetime.date.today().year - \
        all_samples_data["date_of_birth"].dt.year

    # Load the phecode cohort into a dataframe
    phecode_df = pd.read_csv(args.phecode_cohort_file)
    
    # Select the concept id of interest
    phecode_df = phecode_df[phecode_df["phecode"].astype(str) == str(args.concept_id)]
    cases = phecode_df[phecode_df["count"].astype(int) >= args.phecode_count]
    
    # Merge cases with sample data
    case_data = all_samples_data[all_samples_data["person_id"].isin(cases["person_id"])]
    case_data["phenotype"] = 1
    case_data = case_data[['person_id', 'phenotype', "age", "sex_at_birth"]]

    # Gather controls
    # For controls we exclude all sample that have any phecode count in the phecode_count file.
    # This means the phecode_count for all control samples and this phenotype is zero.
    # This is consistent with how PheTK defines controls.
    excluded_from_controls = phecode_df["person_id"]
    controls = all_samples_data[~all_samples_data["person_id"].isin(excluded_from_controls)]
    controls["phenotype"] = 0
    control_data = controls[['person_id', 'phenotype', "age", "sex_at_birth"]]
    
    data = pd.concat([case_data, control_data])
    if args.verbose:
        print("Num cases: ", len(case_data))
        print("Num controls: ", len(control_data))
        print("Num samples: ", len(samples))
        print("Num samples with phecode data: ", len(phecode_df))

    # Rename and infer columns.
    data["sex_at_birth_Male"] = data["sex_at_birth"].apply(lambda x: 1 if x == "Male" else 0)
    
    filename = os.path.join(args.outdir, args.phenotype + "_phenocovar.csv")
    data[['person_id', 'phenotype',
          "age",
          "sex_at_birth_Male"]].to_csv(
        filename,
        index=False, sep=",")

def parse_lab_measurements(samples, demog, args):
    # Process lab measurements
    ptdata = SQLToDF(aou_queries.ConstructTraitSQL(args.concept_id, args.ppi))
    if args.units is None:
        print("Error: --units is required.")
        exit(1)
    data = pd.merge(ptdata, demog, on="person_id", how="inner")
    data["sex_at_birth_Male"] = data["sex_at_birth"].apply(lambda x: 1 if x == "Male" else 0)
    MSG("After merge, have %s data points"%data.shape[0])


    # Filtering
    data.dropna(axis=0, subset=['value_as_number'],inplace=True)
    MSG("After filter NA, have %s data points"%data.shape[0])
    print("units:", ptdata["unit_concept_name"].unique())
    MSG("  Allowable units: %s"%str(aou_queries.GetUnits(args.units)))
    MSG("  Unique units observed: %s"%(str(set(data["unit_concept_name"]))))
    data = data[data["unit_concept_name"].isin(aou_queries.GetUnits(args.units))]
    MSG("After filter units, have %s data points"%data.shape[0])
    if args.range is not None:
        minval, maxval = aou_queries.GetPhenotypeRange(args.range)
        data = data[(data["value_as_number"]>minval) & (data["value_as_number"]<maxval)]
        MSG("After filter range, have %s data points"%data.shape[0])

    # Filter outlier values based on number of SDs
    num_sds = getattr(args, "outlier_sd", None)
    if num_sds is not  None:
        avg = data["value_as_number"].mean()
        sd = data["value_as_number"].std()
        minval = avg - num_sds * sd
        maxval = avg + num_sds * sd
        data = data[(data["value_as_number"] >= minval) & (data["value_as_number"] <= maxval)]
        MSG("After outlier filtering, have %s data points"%data.shape[0])

    # Determine a single representative value per person
    data['Year'] = data['measurement_datetime'].dt.year
    median_per_year = data.groupby(['person_id','Year']).agg(median_year=('value_as_number', my_median)).reset_index()
    median_of_medians = median_per_year.groupby(['person_id']).agg(median_median=('median_year', my_median)).reset_index()
    median_of_medians.rename({"median_median": "value_as_number"}, inplace=True, axis=1)

    # Merge back to whole dataframe to only keep median of median value per person
    filtered = pd.merge(data, median_of_medians, on=["person_id", "value_as_number"])
    MSG("After merge median medians, have %s data points"%filtered.shape[0])

    # De-duplicate to keep one entry per person
    filtered = filtered.sort_values("measurement_datetime").drop_duplicates(subset=["person_id"], keep="last")
    MSG("After dedup, have %s data points"%filtered.shape[0])

    # Output histogram of phenotype values
    plt.hist(filtered["value_as_number"])
    plt.savefig(args.phenotype+"_histogram.png", dpi=300)

    # Record age info
    filtered["age"] = filtered['measurement_datetime'].dt.year - \
        filtered["date_of_birth"].dt.year

    # Optionally add any trait specific drugexposure covariates
    covar_cols = []
    if args.drugexposure_covariate_concept_ids is not None:
        covar_concepts = args.drugexposure_covariate_concept_ids.strip().split(",")
        for ci in covar_concepts:
            concept_id, concept_name = ci.split(":")
            drugdata = GetDrugData(concept_id)
            drugdata = pd.merge(drugdata, filtered[["person_id", "measurement_datetime"]], on="person_id")
            drugdata[concept_name] = drugdata.apply(lambda x: \
                OverlapDrugMeasurement(x["measurement_datetime"], x["start"], x["end"]), 1)
            filtered = pd.merge(filtered, drugdata[["person_id", concept_name]], \
                on="person_id", how="left").fillna(value={concept_name: 0})
            MSG("After add %s, filtered has %s data points"%(concept_name, filtered.shape[0]))
            covar_cols.append(concept_name)

    # Output final phenotype value
    MSG("Final file has %s data points"%filtered.shape[0])
    filtered.rename({"value_as_number": "phenotype"}, inplace=True, axis=1)
    filename = os.path.join(args.outdir, args.phenotype + "_phenocovar.csv")
    filtered[["person_id", "phenotype", "age", "sex_at_birth_Male"]+covar_cols].to_csv(filename, index=False)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotype ID", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs,sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--concept-id", help="Concept ID or comma separated IDs for phenotype.",
                        type=str)
    parser.add_argument("--drugexposure-covariate-concept-ids",
                        help="Comma-separated list of conceptid:conceptname to use as drug exposure covariates", type=str)
    parser.add_argument("--units", help="Comma-separated list of acceptable units. Accepted shorthands: blood", type=str)
    parser.add_argument("--range", help="min, max acceptable phenotype values", type=str)
    parser.add_argument("--ppi", help="Whether or not the phenotype is from the PPI measurements " + \
                                      "as opposed to LOINC. Only physical measurements (e.g. height) " + \
                                      "are in PPI.",
                                      action="store_true", default=False)
    parser.add_argument("--phecode-cohort-file", help="File with phenocodes and counts for each sample." + \
                            "If provided, no queries will be made for phenotype extraction and only the phecode " + \
                            "file is used. --phecode-count should be set to use this file.")
    parser.add_argument("--phecode-count", help="Used only in the case that --phecode-cohort-file is present." + \
                            "Minimum number of phecode counts present for each sample to be considered a case.",
                        type=int)
    parser.add_argument("--verbose", help="Adjust verbosity.",
                                      action="store_true", default=False)
    parser.add_argument("--outlier-sd", help="filter samples with phenotype values exceeding this number of SDs",
                                        type=int, required=False)
    parser.add_argument("--outdir", help="Path to where the pheno_covar file is saved.", type=str)

    args = parser.parse_args()
    
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    if (args.phecode_count is not None and args.phecode_cohort_file is None) or \
        (args.phecode_count is None and args.phecode_cohort_file is not None):
        print("Error: Both --phecode-cohort-file and --phecode-count should be provided.")
        exit(1)
        
    MSG("Processing %s"%args.phenotype)
    # Set up dataframes
    demog = SQLToDF(aou_queries.demographics_sql)

    # Restrict to samples we want to keep
    sampfile = args.samples
    if sampfile.startswith("gs://"):
        sampfile = sampfile.split("/")[-1]
        if not os.path.isfile(sampfile):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(args.samples))
    samples = pd.read_csv(sampfile)

    if args.concept_id and args.phecode_cohort_file is None:
        parse_lab_measurements(
                    samples=samples,
                    demog=demog,
                    args=args,
                    )
    elif args.concept_id and args.phecode_cohort_file:
        parse_phecode_cohort_file(
                    samples=samples,
                    demog=demog,
                    args=args,
                    )

if __name__ == "__main__":
    main()

