#!/bin/bash

# ACAN VNTR chr15:88855424-88857434
# 10Mbp region is chr15:83855424-93857434
# 100Mbp region is chr15:1-100000000
# 1Mbp region is chr15:88355424-89357434
#region="chr15:83855424-93857434"
#region="chr15:1-100000000"
#phenotypes="height"

# INS VNTR
# 1mbp region
#region="chr11:1661570-2661976"
# 10mbp region
#region="chr11:1-12161976"
# whole chromosome
#region="chr11:1-135127769"
# Type 2 Diabetes, 2 phenonames
#phenotypes="EM_202.2 250.2"

# Stroke, CAD and IL1RN VNTR
# IL1RN VNTR chr2:113130529-113130872
#region="chr2:1-244000000" # Whole chromosome
#region="chr2:112630529-113630872" # 1mb region
#region="chr2:108130529-118130872" # 10mb region
#phenotypes="CV_431 CV_431.1 CV_431.11 CV_431.12"


# mean corpuscular hemoglobin concentration and CUL4A VNTR
#CUL4A VNTR chr13:113231980-113232318
#region="chr13:112931980-113532318"
#region="chr13:1-114000000"
#phenotypes="3009744 GE_970"


# atrial fibrillation and NACA VNTR
# NACA VNTR chr12:56717496-56718739
#region="chr12:51717496-61718739" #10mbp region
#phenotypes="CV_416.2 CV_416.21 427.2 427.21 CV_416.211 CV_416.213 CV_416.212 836828 CV_416.214"
#phenotypes="CV_416.2 CV_416.21 427.2 427.21 CV_416.211" # Phenotypes with at least 2000 cases


# Glaucoma and TMCO1 VNTR
# TMCO1 VNTR chr1:165761973-165762288
#region="chr1:160761973-170762288" # 10mbp region
#phenotypes="S01E S01FB S01EA 365 SO_375.1 SO_375.11" # phenotypes with at least 2000 cases, including drugs
#phenotypes="365 SO_375.1 SO_375.11" # phenotypes with at least 2000 cases, excluding drugs
	#--local \

# Colorectal cancer and colon pylops and EIF3H VNTR
# EIF3H VNTR chr8:116622815-116622935
#region="chr8:116122815-117122935" # 10mbp region
#phenotypes="153 CA_101.41" # colorectal cancer 
#phenotypes="153.2" # colon cancer 683 cases
#phenotypes="208 CA_136.4 CA_136.41" # Benight neoplasm of colon
#phenotypes="836844" # colon polyps (PFHH)
#phenotypes="CA_101.411" # Malignant neoplasm of colon
#phenotypes="153 CA_101.41 153.2 208 CA_136.4 CA_136.41 836844 CA_101.411"


# ADHD, PK, SCZ and BP with SLC6A3 VNTR
# SLC6A3 VNTR chr5:1393581-1393986
#region="chr5:1-6393986"
#phenotypes="313.1" # ADHD
#NS_324.1 NS_324.11332 for parkinson's, all below 1000 cases
# 295 295.1 MB_287.1 Schizophrenia
# MB_286.1 296.1 MB_286.11 MB_286.12 836803 bipolar disorder
#phenotypes="313.1 NS_324.1 NS_324.11332 295 295.1 MB_287.1 MB_286.1 296.1 MB_286.11 MB_286.12 836803"

# ADHD, SCZ and OCD with SLC6A4 VNTR
# SLC6A4 VNTR chr17:30221385-30221590
#region="chr17:25221385-35221590"
#phenotypes="313.1" # ADHD
# 295 295.1 MB_287.1 Schizophrenia
# MB_289 300.3 ocd
#phenotypes="313.1 295 295.1 MB_287.1 MB_289 300.3"

# ADHD, OCD with DRD4 VNTR
# DRD4 VNTR chr11:639988-640180
#region="chr11:1-5640180"
#phenotypes="313.1 MB_289 300.3"


# Polyarthritis and IL4 VNTR
# IL4 VNTR chr5:132680584-132680794
#region="chr5:127680584-137680794"
#phenotypes="MS_708 740.1 MS_705 MS_708.1 MS_705.1 714 714.1 836818 MS_708.11" # All phenotypes with at least 1000 cases


# age of onset BD and PER3 VNTR
# PER3 VNTR chr1:7829873-7830144
#region="chr1:2829873-12830144"
#phenotypes="MB_286.1 296.1 MB_286.11 MB_286.12 836803"

# ALS and C9orf72 VNTR
# C9orf72  VNTR chr9:27573484-27573555
#region="chr9:22573484-32573555" #10mbp region
#phenotypes="NS_328 290" # Dementias and cerebral degeneration
#phenotypes="NS_328.1 290.1" # Dementias < 900 individuals
#phenotypes="N06D" # Anti dimentia drugs
#phenotypes="NS_328 290 NS_328.1 290.1 N06D"


# Dystonia and PK and TAF1 VNTR
# TAF1 VNTR chrX:71453055-71453129
#region="chrX:66453055-76453129" #10mbp
#phenotypes="NS_324.31" # dystonia n_cases=869
#phenotypes="NS_324.1 NS_324.11" # parkinson's disease n_cases=751 and 700
#phenotypes="332" # parkinson's disease n_cases=300
#phenotypes="N04" # Anti parkinson's drugs
#phenotypes="N04 NS_324.31 NS_324.1 NS_324.11 332"

# Epilepsy and CSTB VNTR
# CSTB VNTR chr21:43776443-43776479
#region="chr21:38776443-48776479" #10mbp
#phenotypes="NS_330 345 345.1 NS_330.1" # epilepsy n=5300 and 4300 and 1500
#phenotypes="NS_330.12 345.12" # partial epilepsy n=800, 400
# NS_330.11 is 200 for partial epilepsy.
# 836852 is 200 for PFHH
#phenotypes="NS_330 345 345.1 NS_330.1 NS_330.12 345.12"


# Dimentia and PRNP VNTR
# PRNP VNTR chr20:4699397-4699493
#region="chr20:1-9699493" #<10mbp
#phenotypes="NS_328 290 NS_328.1 290.1 N06D"

# motor neuropathy and VWA1
# VWA1 VNTR chr1:1435799-1435818
region="chr1:1-6435818" # 5mbp region
#phenotypes="GE_972.1" # Hereditary neuropathy n=3072
#phenotypes="357" # Inflammatory and toxic neuropathy n=1700
# 250.6 #Polyneuropathy in diabetes n=1568
# 356 Hereditary and idiopathic peripheral neuropathy n_cases=1319
# SO_380.1 Optic neuropathy n=1028
# 836811 PFHH
# NS_337.3 Toxic polyneuropathy n=500
# 377.3 Optic neuritis/neuropathy n=400
# NS_343.1 Peripheral autonomic neuropathy n=265
# NS_337.6 Idiopathic progressive neuropathy n=229
phenotypes="GE_972.1 357 250.6 356 SO_380.1 836811 NS_337.3 377.3"

for phenotype in $phenotypes; do
  python query_allxall_associations.py \
	--type ACAF \
        --pop META \
	--phenoname $phenotype \
	--region $region
done
