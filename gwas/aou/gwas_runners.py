"""
Classes for performing GWAS
"""

import numpy as np

MT_WGS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt' 
SMALLNUM = 10e-400

class HailRunner:
    import hail as hl # only import hail if we have to
    def __init__(self, ptcovar, sample_call_rate = 0.9, variant_call_rate = 0.9, MAF = 0.01, HWE = 1e-15, GQ = 20,region=None, covars=None):
        self.ptcovar = ptcovar
        self.sample_call_rate = sample_call_rate
        self.variant_call_rate = variant_call_rate
        self.MAF = MAF
        self.HWE = HWE
        self.GQ = GQ
        self.region = region
        self.covars = covars  
        self.gwas = None
        self.method = "hail"
        self.setup()

    def setup(self):
    	# Set up hail
        self.hl.init(default_reference = "GRCh38")

        # Load genotypes
        mt = self.hl.read_matrix_table(MT_WGS_PATH)
        if self.region is not None:
            mt = self.hl.filter_intervals(mt, [self.hl.parse_locus_interval(self.region,)])

        # Load phenotype and covariates
        ptcovar = self.hl.Table.from_pandas(self.ptcovar, key="person_id")
        data = mt.annotate_cols(ptcovar = ptcovar[mt.s])

        # Genotype QC
        data = data.annotate_entries(FT = self.hl.coalesce(data.FT,'PASS'))
        data = data.filter_entries(data.FT =='PASS')
        data = data.filter_entries(data.GQ >= self.GQ)  #20

        # Locus and Sample QC
        data = self.hl.variant_qc(data)
        data = self.hl.sample_qc(data)
        data = data.filter_cols(data.sample_qc.call_rate >= self.sample_call_rate, keep = True) #0.9
        data = data.filter_rows(data.variant_qc.call_rate >= self.variant_call_rate, keep = True) #0.9
        data = data.filter_rows(self.hl.min(data.variant_qc.AF) > self.MAF, keep = True) #0.01
        data = data.filter_rows(data.variant_qc.p_value_hwe > self.HWE, keep = True) # 1e-15


         # Keep track of data
        self.data = data

    def RunGWAS(self):
        linear_r = self.hl.linear_regression_rows(
            y= self.data.ptcovar.phenotype,
            x= self.data.GT.n_alt_alleles(),
            covariates = [1.0] + [self.data.ptcovar[item] \
            	for item in self.covars]
        )
        gwas = linear_r.annotate(p_value_str= self.hl.str(linear_r.p_value)).to_pandas()
        gwas["chrom"] = gwas["locus"].apply(lambda x: str(x).split(":")[0])
        gwas["pos"] = gwas["locus"].apply(lambda x: int(str(x).split(":")[1]))
        gwas["p_value"] = gwas.apply(lambda x: float(x["p_value_str"])+SMALLNUM, 1)
        gwas["-log10pvalue"] = -np.log10(gwas["p_value"])
        self.gwas = gwas