#!/bin/bash

# The concept id for colon polyp is in SNOMED vocab
# So the query is likely different from LOINC/PPI.
# I'll have to change the query.

python aou_phenotype_preprocessing.py \
	--phenotype diabetes \
   	--concept-id 201820 \
   	--units "tmp" \
    --snomed \
   	--range 0,1
#SNOMEDCode: 73211009
#OMOP Concept Id: 201820
exit 0
python aou_phenotype_preprocessing.py \
	--phenotype glucose \
   	--concept-id 40795740 \
   	--units "milligram per deciliter" \
   	--range 0,300
exit 0
python aou_phenotype_preprocessing.py \
	--phenotype hemoglobin_a1c \
   	--concept-id 3004410 \
   	--units "percent" \
   	--range 0,100
exit 0
python aou_phenotype_preprocessing.py \
	--phenotype colon_polyp \
   	--concept-id 4285898 \
   	--units "tmp" \
   	--range 0,1
exit 0
python aou_phenotype_preprocessing.py \
	--phenotype height_ppi \
   	--concept-id 903133 \
	--ppi \
   	--units centimeter \
   	--range 0,250 \

