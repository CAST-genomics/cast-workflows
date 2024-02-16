# Phenotype preprocessing for AoU

This script takes as input a phenotype ID and runs preprocessing to create a file with phenotype values and phenotype-specific covariates ready to be used for GWAS or other applications. 

**Note, this must be run on the AoU workbench.** You can use a minimal analysis environment (e.g. 4 CPU, 15GB RAM).

Usage example:

```
./aou_phenotype_preprocessing.py \
   --phenotype ALT \
   --concept-id 37047736 \
   --units blood \
   --range 0,250 \
   [--drugexposure-covariate-concept-ids 21601855:statin] \
   [--sample samplefile.txt] \
   [--ppi]
```

Required arguments: 

* `--phenotype <STR>`: phenotype ID to use to name the output file
* `--concept-id <INT>`: must be obtained from the AoU workbench data browser (TODO: add instructions)
* `--units <STR>`: comma-separated list of units to allow. Shortcut `--units blood` allows units: "IU/L", "No matching concept", "international unit per liter", "no value", "unit per liter". TODO: Might change this based on units for other phenotypes
* `--range <INT,INT>`: Minimum and maximum accepted values for the phenotype

Optional arguments:
* `--sample <FILE>`: file with list of samples to restrict to (one ID per line). Typically this would a list of samples that passed upstream sample-level QC info.
* `--drugexposure-covariate-concept-ids <STR>`: List of conceptid:conceptname to use as drug exposure covariates.
* `--ppi`: Indicate if the phenotype concept id is from the PPI or LOINC (default) tables. PPI are general physical measurements (e.g. hight) and include very limited concept ids.

The output file will is a csv file named `${phenotype}_phenocovar.csv` with columns: "person_id", "phenotype", "age", plus an additional column for each drug exposure named by the "conceptname" provided for that drug. e.g. the example above will have columns "person_id", "phenotype", "age", "statin".
