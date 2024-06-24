#!/bin/bash


python aou_phenotype_preprocessing.py \
	--phenotype height_ppi \
   	--concept-id 903133 \
	--ppi \
   	--units centimeter \
   	--range 0,250 \

# The concept id for colon polyp is in SNOMED vocab
# So the query is likely different from LOINC/PPI.
# I'll have to change the query.
python aou_phenotype_preprocessing.py \
	--phenotype colon_polyp \
   	--concept-id 4285898 \
   	--units "a" \
   	--range 0,1 \
