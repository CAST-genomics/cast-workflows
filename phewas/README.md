# Running Phewas

This directory helps run Phewas on a single target locus and all phenotypes given either the genotyped or imputed calls. For more information on Phewas see [PheTK toolkit](https://github.com/nhgritctran/PheTK).

## Quick start

The first time `run_phewas.py` is called with a new output directory, it will take longer while as it has to query and write the `my_phecode_counts.csv` file. This file stores the number of occurrence of each phecode for each sample. Once created, future runs will automatically read from the stored file. This file can also be passed to the `aou_phenotype_preprocessing.py` script for generating pheno-covar file for binary phenotypes.

To run Phewas on tandem repeats, the genotyped (or imputed) vcf file should be provided with the `--tr-vcf` flag. Here is an example:

```
python run_phewas.py --locus TMCO1 \
                     --chrom chr1 \
                     --start 165761972 \
                     --n-threads 10 \
                     --tr-vcf $tr_vcf \
                     --outdir outputs
```
This will run phewas on the login node (cloud environment). The cromshell workflow is under construction.

** Note: a larger than default memory might be needed for phewas. If your job is killed, try increasing the memory (RAM) on the login node.**

The output is a csv file `phewas_results_${locus}_${start}.csv` with association testing results for each phenotype and a manhattan plot for the phewas.

To save running time, once each intermediate file is generated for a particular locus (e.g. genotypes csv file, phewas file, etc.), it will be loaded for all the future runs with the same locus and same output directory and not generated from scratch again.
