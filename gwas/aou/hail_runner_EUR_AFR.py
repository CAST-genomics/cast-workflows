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

        # Restrict MT to samples in either group (for MAF/HW filtering only)
        ref_mt = mt.filter_cols(hl.is_defined(eur_tbl[mt.s]) | hl.is_defined(afr_tbl[mt.s]))

        # filter multiallelics
        ref_mt = ref_mt.filter_rows(hl.len(ref_mt.alleles) == 2)
        
        # Run variant_qc separately for each group
        ref_mt = ref_mt.annotate_cols(
            eur_cohort = hl.is_defined(eur_tbl[ref_mt.s]),
            afr_cohort = hl.is_defined(afr_tbl[ref_mt.s])
        )

        eur_qc = hl.variant_qc(ref_mt.filter_cols(ref_mt.eur_cohort))
        afr_qc = hl.variant_qc(ref_mt.filter_cols(ref_mt.afr_cohort))

        eur_rows = eur_qc.rows()
        afr_rows = afr_qc.rows()


        # Annotate per-group QC back to the main MT
        mt = mt.annotate_rows(
            eur_AF = eur_rows[mt.row_key].variant_qc.AF,
            eur_HWE = eur_rows[mt.row_key].variant_qc.p_value_hwe,

            afr_AF = afr_rows[mt.row_key].variant_qc.AF,
            afr_HWE = afr_rows[mt.row_key].variant_qc.p_value_hwe,
        )

        # Variant filter: keep if either group passes all thresholds
        mt = mt.filter_rows(
            (
                (hl.min(mt.eur_AF) >= self.MAF) &
                (mt.eur_HWE > self.HWE)
            )
            |
            (
                (hl.min(mt.afr_AF) >= self.MAF) &
                (mt.afr_HWE > self.HWE)
            )
        )
        
        # Load phenotype and covariates
        ptcovar = hl.Table.from_pandas(self.ptcovar, key="person_id")
        mt = mt.annotate_cols(ptcovar = ptcovar[mt.s])

        # now filter samples to given cohort
        ids = pd.DataFrame(self.ptcovar['person_id'])
        sample_tbl = hl.Table.from_pandas(ids, key="person_id")
        mt = mt.filter_cols(hl.is_defined(sample_tbl[mt.s]))

        # Genotype QC
        mt = mt.annotate_entries(FT = hl.coalesce(mt.FT,'PASS'))
        mt = mt.filter_entries(mt.FT =='PASS')
        mt = mt.filter_entries(mt.GQ >= self.GQ) #20

        # sample QC
        mt = hl.sample_qc(mt)
        mt = mt.filter_cols(mt.sample_qc.call_rate >= self.sample_call_rate, keep = True) #0.9
        # remaining variant QC
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.call_rate >= self.variant_call_rate, keep = True) #0.9

         # Keep track of data
        self.data = mt

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
