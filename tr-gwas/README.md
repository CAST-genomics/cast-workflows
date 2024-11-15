# Running TR-GWAS on the All of Us workbench

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

## Run a targeted phenotype and cohort

Example test to run tr-gwas :

```
./tr_gwas_aou.py \
--name test \
--phenotype platelet_count \
--cohort EUR,AFR 

```

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

## Run a full job on all samples and all phenotype 

```
./tr_gwas_aou.py \
--name test_all


```

## Output files 
```
* plink friendly phenotype files
<phenotype>_covar_combined.txt
<phenotype>_pheno_plink.txt

* TR-gwas results 
<phenotype>_<cohort>_gwas.tab
```