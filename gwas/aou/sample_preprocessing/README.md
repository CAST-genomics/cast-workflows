# Sample preprocessing for AoU

This directory contains scripts for performing sample qc. The list of samples that pass qc checks can be found at `${WORKSPACE_BUCKET}/samples/passing_samples_v7.1.csv`. 

The sample qc scripts in this directory should be executed in the following order: 
1. `compute_ratios.py` : obtains #singletons insertions/deletions ratio, heterozygous/homozygous ratio, transition/transversion ratio per sample
2. `merge_batches.py`: the previous script runs in batches, this script combines those for further processing
3. `finish_sampleqc.py`: this uses the values computed in the previous script to perform sample qc and outputs a two-column list (`person_id` and `sex_at_birth_Male`) of samples that pass qc checks

### Running `compute_ratios.py`

This script must be run on the AoU workbench and requires a dataproc cluster for hail/spark. Use the command `nohup python3 compute_ratios.py > output.log 2>&1 &`.

The image belows provides the parameters used to set up the cluster.

![unnamed](https://github.com/CAST-genomics/cast-workflows/assets/16807372/52801a8d-d28b-4123-8011-5161d2a2a950)

With the parameters above, batches of 50k samples are recommended. As there are currently 250k samples available in AoU, the default number of batches is 5.
The `compute_ratios.py` script should, but does not, run through all the batches. The cluster runs out of memory after each batch. The solution, though not ideal, is
to do the following:
- Run one batch. Soon after the next batch starts, it will likely run out of memory.
- Delete and recreate the compute environment.
- In `compute_ratios.py`, change the `for` loop from `for i in range(num_batches):` to `for i in range(X, num_batches):` where `X` is the batch you wish to run. For example, if you have finished two batches and want to run the third batch, `X` would be 2 (following 0 indexing)

Better solutions may be implemented in the future, but as we expect sample preprocessing to be an infrequent step (only performed when new samples are available), we defer this to later.

### Sample qc

After `compute_ratios.py` has finished all batches, run `merge_batches.py` to combine the batches into one file. Then run `finish_sampleqc.py` to perform sample qc and obtain a list of samples that pass qc checks.

The sample preprocessing is performed in `finish_sampleqc.py`. This script uses the ratios to filter samples. 

Related samples are removed as described in https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized

Lastly, this script performs sex-related preprocessing.

# Constructing cohorts

The `construct_cohorts.py` script uses the final output of sample qc above to build ancestry cohorts. For now, it constructs only two cohorts:
1. A European cohort consisting of samples where predicted ancestry is EUR and self-reported race is white
2. An African cohort consiting of samples where predicted ancestry is AFR and self-reported race is black and African American

Each cohort file consists of the same two columns as before: `person_id` and `sex_at_birth_Male`. The cohort files are saved to `${WORKSPACE_BUCKET}/samples/`. Additional cohorts may be added in the future.
