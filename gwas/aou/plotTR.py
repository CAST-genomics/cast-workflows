#!/usr/bin/env python3

"""
Make a plot summarizing TR/trait association and allele freqs

Example:
./plotTR.py --phenotype platelet_count --tr-vcf CBL_test.filtered.sorted.vcf.gz --region chr11:119206290-119206323 --outprefix CBL_platelet_count
"""

# Allow making plots even with no x-forward
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse
import numpy as np
import os
import pandas as pd
import scipy.stats
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.weightstats

from utils import MSG, ERROR

SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)


#use jonanthan's to add CI, returns (mean prob, lower ci bound, upper ci bound)
def weighted_binom_conf(data,weights, successes, confidence):
    assert data[weights].shape == data[successes].shape
    assert len(weights.shape) == 1
    t = np.sum(weights)
    phat = np.dot(weights, successes)/t
    z = scipy.stats.norm.ppf(1 - confidence/2)
    c = z * np.sqrt(np.dot(weights, weights))

    divisor = 2 + 2*c**2/t**2
    center = (2*phat + c**2/t**2)/divisor
    interval_size = c/t*np.sqrt(4*phat*(1-phat) + c**2/t**2)/divisor

    return (phat, center - interval_size, center + interval_size)




def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes.", type=str, required=True)
    parser.add_argument("--region", help="chr:start-end of target TR", type=str, required=True)
    parser.add_argument("--outprefix", help="Output prefix", type=str, required=True)
    parser.add_argument("--min-samples-per-dosage", help="Filter dosages with fewer samples", type=int, default=100)
    args = parser.parse_args()

    # Load phenotype data
    if args.phenotype.endswith(".csv"):
        ptcovar_path = args.phenotype
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)
    data = pd.read_csv(ptcovar_path)
    data["person_id"] = data["person_id"].apply(str)

    # Load TR genotypes
    invcf = utils.LoadSingleReader(args.tr_vcf, checkgz=True)
    samples = invcf.samples
    region = invcf(args.region)
    nrecords = 0
    for record in region:
        trrecord = trh.HarmonizeRecord(trh.VcfTypes["hipstr"], record)
        afreqs = trrecord.GetAlleleFreqs()
        genotypes = trrecord.GetLengthGenotypes()
        dosages = [sum(item)/2 for item in genotypes]
        trdf = pd.DataFrame({"person_id": samples, "tr_dosage": dosages})
        nrecords += 1
    if nrecords == 0:
        ERROR("No matching TR records found")
    if nrecords > 1:
        ERROR("Multiple matching TR records found")

    # Merge phenotype and TR dosages
    df = pd.merge(data, trdf, on=["person_id"])
    pltdata = df.groupby("tr_dosage", as_index=False).agg(phenotype_mean=("phenotype", np.mean), n=("phenotype", len))
    pltdata = pltdata[pltdata["n"]>args.min_samples_per_dosage]

    # compute CI
    df['len_'],df['95_CI_lower'],df['95_CI_upper'] = df.apply(lambda x:  statsmodels.stats.proportion.proportion_confint(method='wilson')(x.tr_dosage,x.phenotype,0.05),axis=1)

    print(df)
    print(pltdata.head())

	# Plot - TODO
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(pltdata["tr_dosage"], pltdata["phenotype_mean"])
    fig.savefig(args.outprefix + ".png")
if __name__ == "__main__":
    main()