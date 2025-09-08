# Running TargetTR on the All of Us workbench

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
cd cast-workflows/long-tr
```

4. Install dependency for configure file
```
python3 -m pip install requests
/opt/conda/bin/python3 ../utilis/configure-cromshell.py

```

## Run a small test job

It is recommended to first run a small test on a couple samples to make sure everything is set up correctly. e.g.:

```
./longtr_launcher_aou.py \
  --tr-bed mr.big.bed \
  --name mytest \
  --batch-size 2 \
  --batch-num 2 \
  --action both
```

This will print out a friendly turtle with the job ID if successful. Use the following to check the status of your job. It will take around 10-20 minutes to run. If successful the status will eventually change to "Succeeded".

```
cromshell status $JOBID
```

If you see "Failed", you can look at the logs to see what happened:

```
cromshell logs -s ALL $JOBID
```

If successful, you can find the output VCF at:
```
${WORKSPACE_BUCKET}/cromwell-execution/target_longTR/${JOBID}/call-sort_index/${JOBNAME}.filtered.sorted.vcf.gz
```

You can copy it to your workspace to check it:
```
gsutil cp ${WORKSPACE_BUCKET}/cromwell-execution/target_longTR/${JOBID}/call-sort_index/${JOBNAME}.filtered.sorted.vcf.gz .
```

## Run a full job on all samples

The above runs the test TR provided in `test.bed` all samples. You will need to change the `--tr-bed` file and `--name` options according to your run.

```
./longtr_launcher_aou.py \
  --tr-bed test.bed \
  --name myrunname \
  --action run-batches 
```

Note for the full job, you might need to add the `-t` option when checking status so the request doesn't time out, e.g.:
```
cromshell -t 200 status $JOBID
```

The job will take 2-3 hours to run on all samples.

## Detailed usage for targetTR_launcher_aou.py

Required options:
* `--name <STRING>`: Name of the run. Used as the prefix to output files.
* `--tr-bed <PATH>`: Path (local file or GCP path) to BED file of TRs to run (HipSTR reference format). 
* `--action <STRING>`: One of: `create-batches`, `run-batches`, or `both`. Typically set to `run-batches`. See "Modifying the batch size" below.

Additional input files:
* `--file-list <PATH>`: GCP path to manifest file describing sequencing data paths. Defaults to the AoU v7 long read manifest file.
* `--genome-id <PATH>`: GCP path to reference genome. Defaults to a path to `Homo_sapiens_assembly38.fasta`.
* `--genome-idx-id <PATH>`: GCP path to the index (`.fai`) file of the reference genome.

Additional run options:
* `--batch-size <INT>`: Number of samples to process in each batch. Default: 300.
* `--batch-num <INT>`: Number of batches to process. Default: -1 (process all batches). This is helpful to set during debugging to consider a small number of batches.
* `--extra-longtr-args <STRING>`: Add more longtr command.
* `--dryrun`: Don't actually submit the cromshell job, just print the command that would have been run.

## Modifying the batch size

**NOTE**: By default we use batches of size 300. We do not recommend changing the batch size. If you do so, you might need to also change settings in the WDL (e.g. memory used by HipSTR).

You can run the following once to precompute batches for a particular batch size:
```
# Note: tr-bed and name are ignored
# when just creating batches
# This takes around 30 minutes
./longtr_launcher_aou.py \
  --tr-bed test.bed \
  --name myrunname \
  --action create-batches --batch-size $batchsize
```

Then, you can use the new batch size in a future run:
```
# Note: tr-bed and name are ignored
# when just creating batches
./longtr_launcher_aou.py \
  --tr-bed test.bed \
  --name myrunname \
  --action run-batches --batch-size $batchsize
```

To run batch creation and the batch jobs at the same time, set `--action both`. This is only recommended for small tests.

## Additional helpful cromshell commands

To check metadata and logs, you can run:
```
cromshell -t 20 metadata $JOBID (increase -t timeout for large batches)
cromshell slim-metadata $JOBID
cromshell logs -s ALL $JOBID

cromshell logs -s ALL -des $JOBID (use -des for descriptive log info)
cromshell logs -s Failed -des $JOBID
```

To list all output files produced by a workflow:
```
cromshell list-outputs $JOBID
```

To get the summarized status of all jobs in the workflow:
```
cromshell count $JOBID
```

# TODO:
* add cost info
