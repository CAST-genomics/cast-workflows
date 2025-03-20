# Imputation UBK
This is a tutorial on how to impute TRs using EnsembleTR reference panel. For details about the EnsembleTR, please check this [repository](https://github.com/gymrek-lab/EnsembleTR).  

# 1.Prerequisites

1. Access to DNA Nexus
* Make sure the [DNA Nexus toolkit](https://documentation.dnanexus.com/downloads) is installed. You can install it with `pip3 install dxpy`
* Run `dx login` to log in with your DNA Nexus credentials.
* Download the jar file for [dxCompiler](https://github.com/dnanexus/dxCompiler/releases). This is a tool that compiles WDLs to convert them into DNA Nexus-compatible workflows.

# 2.Imputation 

First get files that required for imputation: 1) SNP VCFs, 2) SNP-TR reference files, 3) Genetic maps.

* Here, we used the [EnsembleTR reference panel version IV](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v4/ensembletr_refpanel_4_readme.txt). 
Use this command to download the reference panel and upload to DNANexus: 'wget -O - https://\<s3-bucket-url\>/\<file-name\> | dx upload - --destination /project-folder/\<file-name\>'
* Downlaod the [genetic maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip) from Beagle website. Note, the SNP-TR reference panel are use "chr${num}", while the genetic map use "${num}". Make sure add a "chr" prefix to chromosome numbers.

## 2.1 Prepare VCF files

Imputation requires VCF files but UK Biobank only have pVCF (for WGS) or bgen files (for array data). Need tofirst convert those files to VCFs. For pVCFs, there are two version DRAGEN or Graphtyper. One post mention DRAGEN may need further quality control but no details protocol available. 

Here we going to generate VCF from pVCF (graphtyper verson) and bgen files (array data) on chromosome 21 and compare the imputation results.

### 2.1.1 Merge pVCF for imputation

The size of pVCF files are very large, so we need to extract useful field first and then merge those for imputation. A workflow hasbeen created based on [vcf_trimmer](https://github.com/drarwood/vcf_trimmer/tree/master?tab=readme-ov-file) to do this.
* The workflow file is at `../wdl/trim_vcf.wdl`. The size of pvcf before trimming is round 10GB and will be ~3 GB after trimming. The instance requires large disk size to keep all processed files as well as concat those files. Currently, use instance with 16GB memory and 160 GB storage for trimming.

To use this workflow, first run [./compile_wdl.sh](./compile_wdl.sh) to compile the workflow.

The `trim_vcf` will remove tags in pvcf files to reduce the file size and then combine the processed vcfs. It process by chromosomes. Examples are below:

Low quality variatns are removed using `bcftools -i INFO/AAscore>=0.5` which assigned by GraphTyper.

```bash
# run a test
./vcfTrimmer_launcher_ukb.py --chrom chr21  --batch-size 2 --batch-num 2

# run all chromosomes
./vcfTrimmer_launcher_ukb.py --chrom chr21  --concurrent
```

### 2.1.2 Convert bgen files to vcf

Bgen files is located at `Bulk/Genotype Results/Genotype calls`. Need to check any of those variants need to be filtered. The [bgen files](https://github.com/dnanexus/UKB_RAP/blob/main/end_to_end_gwas_phewas/bgens_qc/bgens_qc.wdl) as well as the resultsing [SNPs need QC](https://github.com/dnanexus/UKB_RAP/blob/main/end_to_end_gwas_phewas/run_array_qc.sh) and it is on [b37 coordinates](https://github.com/dnanexus/UKB_RAP/tree/main/end_to_end_gwas_phewas/liftover_plink_beds_tmp). 


### 2.1.3 Split VCF files by samples

To reduce the memory, the VCF files combined by chromosomes are furhter split by samples. Each batch contains 1k samples.
