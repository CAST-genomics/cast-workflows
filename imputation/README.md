# Running imputation on the All of Us workbench

The goal of imputation is to perform phasing and imputation on a cohort of samples. The workflow performs the following steps:

1. **Phase and run imputation with beagle.** Beagle performs multi-sample calling in batches. This is because using beagle to phase and imputate a huge number of samples (e.g. 200,000) is too memory intensive. After imputate each batch, it uses mergeSTR and bcftools to combine the output VCF files.
2. **Merge imputed TR and SNP seperately with MergeSTR and bcftools.** 
3. **Sort, bgzip and tabix index the output VCF file.**


**For most of our use cases, users do not interact with the WDL described here directly, but rather call launcher scripts ** which handle setting up inputs and calling the WDL on particular cloud systems. Sections below give additional WDL details which can be helpful for development/testing or debugging.


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
git checkout batch_imputation2
cd cast-workflows/imputation
./configure-cromshell.py
```

## Run a small test job

It is recommended to first run a small test on a couple samples to make sure everything is set up correctly. e.g.:

```

./batch_imputation_aou.py \
--name batch2_test \
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.bref3 \
--mem 40 \
--batch-num 2 \
--ref $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.vcf.gz \
--region chr21:5101889-5151889 \
--subset-region

```

This will print out a friendly turtle with the job ID if successful. Use the following to check the status of your job. It will take around 10-20 minutes to run. If successful the status will eventually change to "Succeeded".

```
cromshell status $JOBID
```

If you see "Failed", you can look at the logs to see what happened:

```
cromshell logs -s ALL $JOBID
```

You can check the output:
```
cromshell list-outputs $JOBID
```

## Run a full job on all samples

```
./batch_imputation_aou.py \
--name batch2_test \
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.bref3 \
--mem 100 \
--ref $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.vcf.gz 

```

## To check on status on your run
```
cromshell logs -s ALL -des $JOBID (use -des for descriptive log info)
cromshell logs -s Failed -des $JOBID
```
