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
        self.exit = False
        self.setup()

    def setup(self):
    	# Set up hail
        #hl.init(default_reference = "GRCh38")

        # Load genotypes
        mt = hl.read_matrix_table(MT_WGS_PATH)
        if self.region is not None:
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(self.region,)])

        # filter samples to ensure maf filtering applied to cohort
        ids = pd.DataFrame(self.ptcovar['person_id'])
        sample_tbl = hl.Table.from_pandas(ids, key="person_id")
        data = mt.filter_cols(hl.is_defined(sample_tbl[mt.s]))

        # filter multiallelics
        data = data.filter_rows(hl.len(data.alleles) == 2)

        # Load phenotype and covariates
        ptcovar = hl.Table.from_pandas(self.ptcovar, key="person_id")
        data = data.annotate_cols(ptcovar = ptcovar[data.s])

        print(f"Rows: {data.count_rows()}, Columns: {data.count_cols()}")
        if data.count_rows() == 0:
            self.exit = True
            return

        # region-specific PCA
        eigenvalues, scores, loadings = hl.hwe_normalized_pca(data.GT, k=10)
        data = data.annotate_cols(pca = scores[data.s])

        for i in range(10):
            data = data.annotate_cols(**{f'PC{i+1}': data.pca.scores[i]})

        # Genotype QC
        data = data.annotate_entries(FT = hl.coalesce(data.FT,'PASS'))
        data = data.filter_entries(data.FT =='PASS')
        data = data.filter_entries(data.GQ >= self.GQ) #20

        # Locus and Sample QC
        data = self.hl.variant_qc(data)
        data = self.hl.sample_qc(data)
        data = data.filter_cols(data.sample_qc.call_rate >= self.sample_call_rate, keep = True) #0.9
        data = data.filter_rows(data.variant_qc.call_rate >= self.variant_call_rate, keep = True) #0.9
        data = data.filter_rows(self.hl.min(data.variant_qc.AF) > self.MAF, keep = True) #0.01
        data = data.filter_rows(data.variant_qc.p_value_hwe > self.HWE, keep = True) # 1e-100

         # Keep track of data
        self.data = data

    def RunGWAS(self):
        if self.exit:
            return
        if self.isbinary:
            reg = hl.logistic_regression_rows(
                test = 'wald',
                y= self.data.ptcovar.phenotype,
                x= self.data.GT.n_alt_alleles(),
                covariates = [1.0] + [self.data.ptcovar[item] \
            	    for item in self.covars]
            )
        else:
            pc_covars = [self.data[f'PC{i+1}'] for i in range(10)]
            reg = hl.linear_regression_rows(
                y= self.data.ptcovar.phenotype,
                x= self.data.GT.n_alt_alleles(),
                covariates = [1.0] + [self.data.ptcovar[item] \
            	    for item in self.covars] + pc_covars
            )
        gwas = reg.annotate(p_value_str= hl.str(reg.p_value)).to_pandas()
        gwas["chrom"] = gwas["locus"].apply(lambda x: str(x).split(":")[0])
        gwas["pos"] = gwas["locus"].apply(lambda x: int(str(x).split(":")[1]))
        gwas["p_value"] = gwas.apply(lambda x: float(x["p_value_str"])+SMALLNUM, 1)
        gwas["-log10pvalue"] = -np.log10(gwas["p_value"])
        self.gwas = gwas
