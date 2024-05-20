import pandas as pd
import os
import subprocess
import sys

bucket = os.getenv('WORKSPACE_BUCKET')

subprocess.run("gsutil cp "+bucket+"/phenotypes/creatinine_phenocovar.csv ./", shell=True, check=True)
creatinine = pd.read_csv("creatinine_phenocovar.csv")

subprocess.run("gsutil cp "+bucket+"/samples/AFR_BLACK.csv ./", shell=True, check=True)
subprocess.run("gsutil cp "+bucket+"/samples/NOT_AFR_BLACK.csv ./", shell=True, check=True)
afr = pd.read_csv("AFR_BLACK.csv")
afr['isblack'] = 1
not_afr = pd.read_csv("NOT_AFR_BLACK.csv")
not_afr['isblack'] = 0

afr_ancestry = pd.concat([afr, not_afr], axis=0).reset_index(drop=True)

allinfo = creatinine.merge(afr_ancestry, on='person_id', how='inner')
allinfo.drop(columns='sex_at_birth_Male_y', inplace=True)
allinfo.rename(columns={'sex_at_birth_Male_x':'sex_at_birth_Male'}, inplace=True)


def egfr_calc(scr, ismale, isblack, age):
    k = 0.9 if ismale else 0.7
    a = -0.411 if ismale else -0.429
    egfr = 141 * (min(scr / k, 1)**a) * (max(scr / k, 1)**-1.209) * (0.993**age)
    if not ismale:
        egfr *= 1.1018
    if isblack:
        egfr *= 1.159
    return egfr

allinfo['egfr'] = allinfo.apply(lambda x: egfr_calc(x.phenotype, x.sex_at_birth_Male, x.isblack, x.age), axis = 1)

afr_egfr = allinfo[allinfo['isblack']==1][['person_id', 'egfr', 'age', 'sex_at_birth_Male']]
afr_egfr.rename(columns={'egfr':'phenotype'}, inplace=True)
not_afr_egfr = allinfo[allinfo['isblack']==0][['person_id', 'egfr', 'age', 'sex_at_birth_Male']]
not_afr_egfr.rename(columns={'egfr':'phenotype'}, inplace=True)

all_egfr = pd.concat([afr_egfr, not_afr_egfr], axis=0).reset_index(drop=True)
all_egfr.to_csv("egfr_ckdepi_phenocovar.csv", index=False) #meaning from the ckd-epi formula bc egfr_phenocovar.csv exists


subprocess.run("gsutil cp egfr_ckdepi_phenocovar.csv "+bucket+"/phenotypes/egfr_ckdepi_phenocovar.csv", shell=True, check=True)
