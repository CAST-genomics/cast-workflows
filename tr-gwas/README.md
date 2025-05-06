# Running TR-GWAS genome-wide on the All of Us workbench

The goal of this workflow is to run TR-gwas genome-wide with plink2 across different phenotypes in different cohorts. 

For most of our use cases, users do not interact with the WDL described here directly, but rather call launcher scripts which handle setting up inputs and calling the WDL using cromshell. Sections below give additional WDL details which can be helpful for development/testing or debugging.


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
cd cast-workflows/tr-gwas
../utils/configure-cromshell.py
```

## Run TR-Gwas on targeted phenotype and cohorts

Example test to run tr-gwas on quantative phenotype :

```
./tr_gwas_aou.py \
--name test_quantative_trait \
--phenotype platelet_count \
--cohort EUR,AFR,AMR 
```

Example test to run tr-gwas on case/control phenotype :

```
./tr_gwas_aou.py \
--name test_binary_trait \
--phenotype platelet_count \
--cohort ALL
--logistic
```

Example test to run tr-gwas on specific ancestry pc :

```
./tr_gwas_aou.py \
--name test_binary_trait \
--phenotype platelet_count \
--cohort EUR \
--ancestry-pc \
--ancestry-pc-path ${filename}

```

Required options:
* `--name <STRING>`: Name of the run. Used as the prefix to output files.
* `--phenotype <STRING>`: Name of the phenotype
* `--cohort <STRING>`: Choose one or more from `EUR`, `AFR`, `AMR`, `NOT_AFR`, `ALL` . If empty, default run all. 

Additional input files:
* `--ancestry-pc-path <STRING>`: Name of the ancestry pc file. Need to upload to google bucket first $WORKSPACE/ancestry_pc/

Additional run options:
* `--logistic`:Run logistic regression. Default is linear regression
* `--ancestry-pc`: Use ancestry specific PC instead of Global PC. Default use global PC
* `--covar-name <STRING>`: Name of covariates to run, default is "age,sex_at_birth_Male,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"



This will print out a friendly turtle with the job ID if successful. Use the following to check the status of your job. It will take around 10-20 minutes to run. If successful the status will eventually change to "Succeeded".

```
cromshell status $JOBID
```

If you see "Failed", you can look at the logs to see what happened:

```
cromshell logs -s ALL $JOBID
cromshell logs -s ALL -des $JOBID (use -des for descriptive log info)
cromshell logs -s Failed -des $JOBID
```

You can check the output:
```
cromshell list-outputs $JOBID
```

## Run a full job on all samples and all phenotypes 

```
./tr_gwas_aou.py \
--name test_all
```

Warning: This script manually filters PC6 for EUR and NOT_AFR cohorts


## Samples preprocessed in cast-workflows/gwas/aou/sample_preprocessing

## Output files 

plink friendly phenotype files
```
<phenotype>_covar_combined.txt
<phenotype>_pheno_plink.txt
```

TR-gwas results 
```
<phenotype>_<cohort>_gwas.tab
```

## Plotting Manhattan plots 

Example to plot manhattan and qq plots:

```
python plot_gwas.py --gwas platelet_count_gwas.tab

## TODO

 (1) looks at pairwise correlation of all covariates + the phenotype and removes highly correlated ones and 
 (2) filtering for any categorical covariates as we do in Margoliash et al. e.g. statin usage and sex are examples of categorial covariates we have used.