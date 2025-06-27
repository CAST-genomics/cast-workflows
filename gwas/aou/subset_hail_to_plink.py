import sys
import os
import hail as hl
hl.init(default_reference = "GRCh38")

bucket = os.getenv('WORKSPACE_BUCKET')
chrom = sys.argv[1]
from_bp = sys.argv[2]
to_bp = sys.argv[3]

outfile = bucket+"/nichole/plink/"+chrom+"_"+from_bp+"_"+to_bp+"_plink"

region = "chr"+chrom+":"+from_bp+"-"+to_bp

mt_wgs_path = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
mt = hl.read_matrix_table(mt_wgs_path)
mt = hl.filter_intervals(mt, [hl.parse_locus_interval(region,)])

# Keep only bi-allelic variants
mt = mt.filter_rows(hl.len(mt.alleles) == 2, keep=True)

# Compute variant QC metrics to get MAF
mt = hl.variant_qc(mt)

# Filter to MAF > 5%
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.05, keep=True)

hl.export_plink(mt, outfile, ind_id = mt.s)