#!/bin/bash

# ACAN VNTR chr15:88855424-88857434
# 10Mbp region is chr15:83855424-93857434
# 100Mbp region is chr15:1-100000000
# 1Mbp region is chr15:88355424-89357434
#region="chr15:83855424-93857434"
#region="chr15:1-100000000"
#phenotype="height"

# INS VNTR
# 1mbp region
#region="chr11:1661570-2661976"
# 10mbp region
#region="chr11:1-12161976"
# whole chromosome
#region="chr11:1-135127769"
# Type 2 Diabetes, 2 phenonames
#for phenotype in"EM_202.2" "250.2"; do

# mean corpuscular hemoglobin concentration and CUL4A VNTR
#CUL4A VNTR chr13:113231980-113232318
#region="chr13:113337980-113272318"
region="chr13:1-114000000"
for phenotype in "3009744" "GE_970"; do

  python query_allxall_associations.py \
	--type ACAF \
        --pop META \
	--phenoname $phenotype \
	--local \
	--region $region
done
exit 0
python query_allxall_associations.py \
	--type ACAF \
        --pop EUR \
	--phenoname $phenotype \
	--region $region
python query_allxall_associations.py \
	--type ACAF \
        --pop AFR \
	--phenoname $phenotype \
	--region $region
