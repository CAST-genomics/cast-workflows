# Pre-compute batches of ACAF VCF files for use in future workflows

This set of workflows precomputes small subsets of the ACAF VCF files. 

The steps below launch the following Workflows with cromshell:

* `wdl/subset_vcf.wdl`: Launch a WDL that creates 1 smaller VCF per sample batch per region
* `wdl/concatenate_batch_vcfs.wdl`: Oragnize all the subset VCF files per batch in a single folder in our workspace for future use

## Subset VCF: Launch jobs

Before running you will need to:
* Run the tasks under "Setup" below (once done these shouldn't need to be rerun)
* Start an All of Us analysis environment and go to the terminal
* Start the All of Us cromshell environment
* From the analysis environment terminal, clone this repository and configure cromshell:

```
git clone https://github.com/CAST-genomics/cast-workflows/
cd cast-workflows/subset_vcf
../utils/configure-cromshell.py
```

To launch a job for one chromosome (Note: these are big jobs. only run one chrom at a time).

```
# Run subset VCF
chrom=XXX # e.g. chrom=11
./subset_vcf_launcher.py --chrom ${chrom} --name chr${chrom} # copy the jobid to $jobid
cromshell alias $jobid subset-vcf-chr${chrom}

# Concatenate the subsets
./concatenate_batches_v2.py subset-vcf-chr${chrom} chr${chrom} # copy the jobid to $concatjobid
cromshell alias $concatjobid concat-vcf-chr${chrom}

# Upload to ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom}
cromshell -mc list-outputs -j concat-vcf-chr${chrom} | \
	python -c "import json, sys; data=json.load(sys.stdin); [sys.stdout.write(item+'\n') for item in data['concatenate_batch_vcfs.vcf_outputs']+data['concatenate_batch_vcfs.vcf_indices']]" | \
	xargs -n1 -I% -P1 sh -c "gsutil mv % ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom}/"
```

Note: large chroms need to be concatenated in batches. This is bc the JSON is too large to run the full chrom. e.g.:

```
./concatenate_batches_v2.py subset-vcf-chr${chrom} chr${chrom} 0 100
./concatenate_batches_v2.py subset-vcf-chr${chrom} chr${chrom} 100 246
```

## Setup

The steps below:

Setup work on the terminal:

* Break up samples into batches of 1000
* Sets up the "groups" file needed for subsetting with bcftools split (`aou_sample_groups.txt`) 
* Sets up 10Mb regions to run on one at a time

### Set up group batches

```
# Set up groups file for batches of 1000
gsutil cp $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample/aou_sample_list.txt .
mkdir batches
split -l 1000 aou_sample_list.txt batches/batch
counter=1
for file in batches/batch*
do
    mv "$file" "batches/batch$counter"
    ((counter++))
done

# Make full group set
for file in batches/batch*
do
	batchname=$(basename $file)
	cat $file | awk -v"prefix=$batchname" '{print $1 "\t-\t"prefix}'
done > aou_sample_groups.txt
gsutil cp aou_sample_groups.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/aou_sample_groups.txt

# Make small group set
cat aou_sample_groups.txt | awk '($3=="batch1" || $3=="batch2")' > aou_sample_groups_test.txt
gsutil cp aou_sample_groups_test.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/aou_sample_groups_test.txt
```

### Set up region batches
```
bedtools makewindows -g hg38.txt -w 10000000 > genome_windows.bed
for chrom in $(seq 1 22)
do
	cat genome_windows.bed | grep -w "chr"${chrom} | awk '{print $1 ":" $2 "-"$3}' > regions_chr${chrom}.txt
	gsutil cp regions_chr${chrom}.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/regions_chr${chrom}.txt
done
```

### Fixing chr9

chr9 last batch had failed so we used the following steps to manually fix:

```
# Fix the batch
bcftools view -Oz -o tmp.vcf.gz -r chr9:130000000-138394717 ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr9.vcf.bgz
tabix -f -p vcf tmp.vcf.gz
bcftools plugin split tmp.vcf.gz -G renamed_sample_groups.txt -Oz -o .

for f in batch*.vcf.gz
do
	tabix -p vcf $f
	gsutil cp $f ${WORKSPACE_BUCKET}/subset_vcf/fix_chr9/
	gsutil cp $f.tbi ${WORKSPACE_BUCKET}/subset_vcf/fix_chr9/
done

# Update the chr9 concat json file
./misc/fix_batches_chr9.py

# Launch
cromshell submit \
	wdl/concatenate_batch_vcfs.wdl \
	chr9concat.fix.aou.json \
	-op chr9concat.aou.options.json
```