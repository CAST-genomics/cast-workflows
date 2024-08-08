#!/bin/bash

# ACAN VNTR chr15:88855424-88857434
# 10Mbp region is chr15:83855424-93857434
# 100Mbp region is chr15:1-100000000
# 1Mbp region is chr15:88355424-89357434
#region="chr15:83855424-93857434"
#region="chr15:1-100000000"

# INS VNTR
# 1mbp region
region="chr11:1661570-2661976"
# 10mbp region
region="chr11:1-12161976"

#CUL4A VNTR chr13:113231980-113232318
region="chr13:1-114000000"

phenotype="height"

# Type 2 Diabetes, 2 phenonames
#phenotype="EM_202.2"
phenotype="250.2"

# mean corpuscular hemoglobin concentration
phenotype="3009744"
python query_allxall_associations.py \
	--type ACAF \
        --pop META \
	--phenoname $phenotype \
	--region $region
phenotype="GE_970"
python query_allxall_associations.py \
	--type ACAF \
        --pop META \
	--phenoname $phenotype \
	--region $region
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
