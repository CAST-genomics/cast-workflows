import pandas as pd
import os
import subprocess
import sys

bucket = os.getenv('WORKSPACE_BUCKET')

subprocess.run("gsutil cp "+bucket+"/tmirmira/data/qc/wgs_ratios_all.csv ./", shell=True, check=True)
all_ratios = pd.read_csv("wgs_ratios_all.csv")

#ratio filtering here
means = all_ratios.mean()
stds = all_ratios.std()
for c in all_ratios.columns:
    if c == 'person_id':
        continue
    avg = means[c]
    sd = stds[c]
    minval = avg - 8*sd
    maxval = avg + 8*sd
    all_ratios = all_ratios[(all_ratios[c] >= minval) & (all_ratios[c] <= maxval)]


#relatedness filtering
related_samples_path = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv"
subprocess.run("gsutil -u $GOOGLE_PROJECT -m cp "+related_samples_path+" ./", shell=True, check=True)
related = pd.read_csv("relatedness_flagged_samples.tsv", sep="\t")
all_ratios = all_ratios[~all_ratios['s'].isin(related['sample_id'].values)]

#v7.1 update as per https://support.researchallofus.org/hc/en-us/articles/25646444720404-Incremental-Data-Release-for-v7-Genotyping-Array-and-Short-Read-Genomic-Data
wgs_to_filter_path = "gs://fc-aou-datasets-controlled/v7/known_issues/research_id_v7_wgs_known_issue_15.tsv"
subprocess.run("gsutil -u $GOOGLE_PROJECT cp "+wgs_to_filter_path+" ./", shell=True, check=True)
wgs_to_filter = pd.read_csv("research_id_v7_wgs_known_issue_15.tsv", sep="\t")
all_ratios = all_ratios[~all_ratios['s'].isin(wgs_to_filter['research_id'].values)]

#sex filtering/processing
# This query represents dataset "All_population" for domain "person" and was generated for All of Us Controlled Tier Dataset v7
dataset_52206452_person_sql = """
    SELECT
        person.person_id,
        person.gender_concept_id,
        p_gender_concept.concept_name as gender,
        person.birth_datetime as date_of_birth,
        person.race_concept_id,
        p_race_concept.concept_name as race,
        person.ethnicity_concept_id,
        p_ethnicity_concept.concept_name as ethnicity,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.person` person
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_gender_concept
            ON person.gender_concept_id = p_gender_concept.concept_id
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_race_concept
            ON person.race_concept_id = p_race_concept.concept_id
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_ethnicity_concept
            ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_sex_at_birth_concept
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
    WHERE
        person.PERSON_ID IN (
            SELECT
                distinct person_id
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person
            WHERE
                cb_search_person.person_id IN (
                    SELECT
                        person_id
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p
                    WHERE
                        has_whole_genome_variant = 1
                )
            )"""

demog = pd.read_gbq(
    dataset_52206452_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

demog = demog[(demog['sex_at_birth']=='Male') | (demog['sex_at_birth']=='Female')]

all_ratios = all_ratios.merge(demog, left_on='s', right_on='person_id', how='inner')
all_ratios = pd.get_dummies(all_ratios, columns=['sex_at_birth'], dtype='int')
all_ratios[['person_id','sex_at_birth_Male']].to_csv("postqc_samples.csv", index=False)
subprocess.run("gsutil cp postqc_samples.csv "+bucket+"/samples/passing_samples_v7.1.csv", shell=True, check=True)
