#!/bin/bash

docker run -v /storage/mgymrek/workspace/hipstr-gcs/hipstr:/ref -v ${PWD}:/data gcr.io/ucsd-medicine-cast/hipstr-gymreklab-gcs HipSTR \
    --bams gs://gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.cram \
    --fasta /ref/Homo_sapiens_assembly38.fasta \
    --regions /data/test.bed \
    --def-stutter-model \
    --min-reads 2 \
    --str-vcf /data/test.vcf.gz
