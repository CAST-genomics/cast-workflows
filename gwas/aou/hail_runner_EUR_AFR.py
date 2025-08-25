"""
Classes for performing GWAS
"""

import hail as hl
import numpy as np
import os
import pandas as pd

MT_WGS_PATH = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
SMALLNUM = 10e-400

class HailRunner:

    import hail as hl
    def __init__(self, ptcovar, isbinary=False, region=None, covars=[], sample_call_rate=None, variant_call_rate=None, MAF=None, HWE=None, GQ=None):

        self.ptcovar = ptcovar
        self.isbinary = isbinary
        self.region = region
        self.covars = covars 
        self.sample_call_rate = sample_call_rate
        self.variant_call_rate = variant_call_rate
        self.MAF = MAF
        self.HWE = HWE
        self.GQ = GQ  
        self.gwas = None
        self.data = None
        self.method = "hail"
        self.setup()

    def setup(self):
    	# Set up hail
        hl.init(default_reference = "GRCh38")

        # Load genotypes
        mt = hl.read_matrix_table(MT_WGS_PATH)
        if self.region is not None:
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(self.region,)])

        os.system("gsutil -u ${GOOGLE_PROJECT} cp ${WORKSPACE_BUCKET}/samples/EUR_WHITE.csv .")
        os.system("gsutil -u ${GOOGLE_PROJECT} cp ${WORKSPACE_BUCKET}/samples/AFR_BLACK.csv .")
        eur_sample_ids = pd.read_csv("EUR_WHITE.csv")['person_id'].astype(str)
        afr_sample_ids = pd.read_csv("AFR_BLACK.csv")['person_id'].astype(str)
        eur_tbl = hl.Table.from_pandas(pd.DataFrame(eur_sample_ids), key="person_id")
        afr_tbl = hl.Table.from_pandas(pd.DataFrame(afr_sample_ids), key="person_id")

        data = mt.filter_cols(hl.is_defined(eur_tbl[mt.s]) | hl.is_defined(afr_tbl[mt.s]))


        # filter multiallelics
        data = data.filter_rows(hl.len(data.alleles) == 2)

        # Load phenotype and covariates
        ptcovar = hl.Table.from_pandas(self.ptcovar, key="person_id")
        data = data.annotate_cols(ptcovar = ptcovar[data.s])


        # Genotype QC
        data = data.annotate_entries(FT = hl.coalesce(data.FT,'PASS'))
        data = data.filter_entries(data.FT =='PASS')
        data = data.filter_entries(data.GQ >= self.GQ) #20

        
        # Run variant_qc separately for each group
        data = data.annotate_cols(
            eur_cohort = hl.is_defined(eur_tbl[data.s]),
            afr_cohort = hl.is_defined(afr_tbl[data.s])
        )

        eur_qc = hl.variant_qc(data.filter_cols(data.eur_cohort))
        afr_qc = hl.variant_qc(data.filter_cols(data.afr_cohort))

        # Annotate per-group QC back to the main MT
        eur_rows = eur_qc.rows()
        afr_rows = afr_qc.rows()

        data = data.annotate_rows(
            eur_AF = eur_rows[data.row_key].variant_qc.AF,
            eur_HWE = eur_rows[data.row_key].variant_qc.p_value_hwe,

            afr_AF = afr_rows[data.row_key].variant_qc.AF,
            afr_HWE = afr_rows[data.row_key].variant_qc.p_value_hwe,
        )

        # Variant filter: keep if either group passes all thresholds
        data = data.filter_rows(
            (
                (hl.min(data.eur_AF) >= self.MAF) &
                (data.eur_HWE > self.HWE)
            )
            |
            (
                (hl.min(data.afr_AF) >= self.MAF) &
                (data.afr_HWE > self.HWE)
            )
        )
        

        # now filter samples to given cohort
        ids = pd.DataFrame(self.ptcovar['person_id'])
        sample_tbl = hl.Table.from_pandas(ids, key="person_id")
        data = data.filter_cols(hl.is_defined(sample_tbl[data.s]))

        # sample QC
        data = hl.sample_qc(data)
        data = data.filter_cols(data.sample_qc.call_rate >= self.sample_call_rate, keep = True) #0.9
        # remaining variant QC
        data = hl.variant_qc(data)
        data = data.filter_rows(data.variant_qc.call_rate >= self.variant_call_rate, keep = True) #0.9

         # Keep track of data
        self.data = data

    def RunGWAS(self):
        if self.isbinary:
            reg = hl.logistic_regression_rows(
                test = 'wald',
                y= self.data.ptcovar.phenotype,
                x= self.data.GT.n_alt_alleles(),
                covariates = [1.0] + [self.data.ptcovar[item] \
            	    for item in self.covars]
            )
        else:
            reg = hl.linear_regression_rows(
                y= self.data.ptcovar.phenotype,
                x= self.data.GT.n_alt_alleles(),
                covariates = [1.0] + [self.data.ptcovar[item] \
            	    for item in self.covars]
            )
        gwas = reg.annotate(p_value_str= hl.str(reg.p_value)).to_pandas()
        gwas["chrom"] = gwas["locus"].apply(lambda x: str(x).split(":")[0])
        gwas["pos"] = gwas["locus"].apply(lambda x: int(str(x).split(":")[1]))
        gwas["p_value"] = gwas.apply(lambda x: float(x["p_value_str"])+SMALLNUM, 1)
        gwas["-log10pvalue"] = -np.log10(gwas["p_value"])
        self.gwas = gwas
