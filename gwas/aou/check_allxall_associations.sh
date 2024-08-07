#!/bin/bash

# ACAN VNTR chr15:88855424-88857434
# 10Mbp region is chr15:83855424-93857434
# 100Mbp region is chr15:1-100000000
# 1Mbp region is chr15:88355424-89357434
region="chr15:83855424-93857434"
region="chr15:1-100000000"
python check_allxall_associations.py \
	--type ACAF \
        --pop META \
	--phenoname height \
	--region $region
