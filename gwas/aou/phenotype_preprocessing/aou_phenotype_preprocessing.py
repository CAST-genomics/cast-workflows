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

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotype ID", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs,sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--concept-id", help="Concept ID or comma separated IDs for phenotype.",
                        type=str, required=True)
    parser.add_argument("--skip-concepts-in-controls", help="Skip controls with this concept ids." + \
                        " Only taken into account for binar (snomed) phenotypes.",
                        type=str)
    parser.add_argument("--drugexposure-covariate-concept-ids",
                        help="Comma-separated list of conceptid:conceptname to use as drug exposure covariates", type=str)
    parser.add_argument("--units", help="Comma-separated list of acceptable units. Accepted shorthands: blood", type=str)
    parser.add_argument("--range", help="min, max acceptable phenotype values", type=str)
    parser.add_argument("--ppi", help="Whether or not the phenotype is from the PPI measurements " + \
                                      "as opposed to LOINC. Only physical measurements (e.g. height) " + \
                                      "are in PPI.",
                                      action="store_true", default=False)
    parser.add_argument("--snomed", help="Whether or not thhe phenotype is from the SNOMED vocab, including diseases.",
                                      action="store_true", default=False)
    parser.add_argument("--verbose", help="Adjust verbosity.",
                                      action="store_true", default=False)
    parser.add_argument("--outlier-sd", help="filter samples with phenotype values exceeding this number of SDs",
                                        type=int, required=False)

    args = parser.parse_args()
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

    concept_ids_int = [int(concept) for concept in args.concept_id.split(",")]
    # Process snomed conditions
    if args.snomed:
        if args.skip_concepts_in_controls:
            all_concepts_for_query = ",".join([args.concept_id, args.skip_concepts_in_controls])
            skip_concepts_in_controls_int = [int(concept) for concept in args.skip_concepts_in_controls(",")]
        # Create the conditions dataframe based on the query
        ptdata = SQLToDF(aou_queries.ConstructSnomedSQL(args.concept_id))
        MSG("Right after query, have %s data points"%ptdata.shape[0])

        # Merge with sample data
        data = pd.merge(ptdata, demog, on="person_id", how="right")
        data = pd.merge(data, samples, on="person_id", how="right")
        
        # Set age and sex columns
        #data["age"] = data.apply(lambda row:
        #        row["condition_start_datetime"].year - row["date_of_birth"].year \
        #        if ((not pd.isnull(row["condition_start_datetime"])) and \
        #            int(row["condition_concept_id"]) == int(args.concept_id)) else\
        #        datetime.date.today().year - row["date_of_birth"].year,
        #    axis=1)
        data["age"] = data.apply(lambda row:
                datetime.date.today().year - row["date_of_birth"].year,
            axis=1)
        MSG("After merge with demog and filter samples, have %s data points"%data.shape[0])
        
        # Compute the binary phenotype column
        # Here we remove samples with concept ids indicated in --skip-concepts-in-controls.
        data = data.fillna(np.nan)
        data["phenotype"] = data["condition_concept_id"].apply(
            lambda x: 1 if not pd.isnull(x) and int(x) in concept_ids_int \
            else 0)
            #else -1 if not pd.isnull(x) and int(x) in skip_concepts_in_controls_int \

        # Merge and groupby to remove duplicate/uninformative rows
        columns = ["person_id", "age", "sex_at_birth_Male", "phenotype"]
        data = data[columns].groupby(["person_id"]).agg(
                        {"phenotype": max,
                        "age": max,
                        "sex_at_birth_Male": \
                            lambda series: \
                            print("Warning: mismatching sex_at_birth in entries of the same person") \
                            if len(set(series)) > 1 else min(series),
                        }
                        ).reset_index()
        if args.verbose:
            print("Number of elements after group_by: ", len(data))
            print("Number of rows with snomed phenotype true: ", len(data[data["phenotype"] == 1]))
            print("Number of rows with snomed phenotype false: ", len(data[data["phenotype"] == 0]))

        # Store the results in the output file
        data[['person_id', 'phenotype',
              "age",
              "sex_at_birth_Male"]].to_csv(
                    args.phenotype+"_phenocovar.csv",
                    index=False, sep=",")
        return
    else:
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
    filtered[["person_id", "phenotype", "age", "sex_at_birth_Male"]+covar_cols].to_csv(args.phenotype+"_phenocovar.csv", index=False)

if __name__ == "__main__":
    main()

