import pandas as pd
import os
import subprocess

bucket = os.getenv('WORKSPACE_BUCKET')

subprocess.run("gsutil cp "+bucket+"/samples/passing_samples_v7.1.csv ./", shell=True, check=True)

passing_samples = pd.read_csv("passing_samples_v7.1.csv")
#columns = person_id, sex_at_birth_Male

ancestry_path = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
subprocess.run("gsutil -u $GOOGLE_PROJECT -m cp "+ancestry_path+" ./", shell=True, check=True)

ancestry = pd.read_csv("ancestry_preds.tsv", sep="\t")
ancestry = ancestry[['research_id', 'ancestry_pred', 'ancestry_pred_other']]
ancestry.rename(columns={'research_id':'person_id'}, inplace=True)

passing_samples = passing_samples.merge(ancestry, on='person_id', how='inner')

dataset_37150131_person_sql = """
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
    dataset_37150131_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

demog = demog[['person_id', 'race']]
passing_samples = passing_samples.merge(demog, on='person_id', how='inner')

eur_white = passing_samples[(passing_samples['ancestry_pred']=='eur') & (passing_samples['race']=='White')]
afr_black = passing_samples[(passing_samples['ancestry_pred']=='afr') & (passing_samples['race']=='Black or African American')]
not_afr_black = passing_samples[~((passing_samples['ancestry_pred']=='afr') | (passing_samples['race']=='Black or African American'))]


eur_white[['person_id', 'sex_at_birth_Male']].to_csv("EUR_WHITE.csv", index=False)
afr_black[['person_id', 'sex_at_birth_Male']].to_csv("AFR_BLACK.csv", index=False)
not_afr_black[['person_id', 'sex_at_birth_Male']].to_csv("NOT_AFR_BLACK.csv", index=False)


subprocess.run("gsutil cp EUR_WHITE.csv "+bucket+"/samples/EUR_WHITE.csv", shell=True, check=True)
subprocess.run("gsutil cp AFR_BLACK.csv "+bucket+"/samples/AFR_BLACK.csv", shell=True, check=True)
subprocess.run("gsutil cp NOT_AFR_BLACK.csv "+bucket+"/samples/NOT_AFR_BLACK.csv", shell=True, check=True)



