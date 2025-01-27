# Extract TR imputation results on the All of Us workbench

To extract specified TR imputation results to perform concordance test between targeTR(Hipstr) and imputation calls

For most of our use cases, users do not interact with the WDL described here directly, but rather call launcher scripts which handle setting up inputs and calling the WDL using cromshell. Sections below give additional WDL details which can be helpful for development/testing or debugging.


The outputs for each chromosome are

${chrom}.vcf.gz
${chrom}.vcf.gz.tbi


## Setup
In all cases you'll have to run the following steps in the AoU workbench before starting a workflow:

1. Start the cromwell environment (pink circle icon on the right).
2. Start a cloud environment (Jupyter icon) with:
    * "General Analysis" environment
    * 2 CPU/7.5GB RAM
    * Compute type: Standard VM
    * Automatically pause after: 30 minutes
    * Storage disk options: standard disk, 120GB
3. Open a terminal (terminal icon) and run the commands below:

```
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/tr-imputation/compare_concordance/
../utils/configure-cromshell.py
```

## Extract ukb finemappped str

It is recommended to first run a small test on a couple samples to make sure everything is set up correctly. e.g.:

```
./run_str_extraction.py \
--name ukb_imputation \
--str ukb_finemapped_hg38_str.txt

```

## Use script to preprocess TGT field and calculate copy num sum for each sample
```
for i in {1..22};do ./preprocess_imputation.sh "ukb_imputation_chr${i}.vcf.gz" "imputation_chr${i}";done 

```

## Use jupyter notebook to process downstream analysis

Imputation_hipstr_comparison_all