#!/bin/bash
VCF=$1
SAMPFILE=$2
OUTPREFIX=$3
GOOGLEPROJECT=$4
# TODO: option to take in a region

# Get column subset to extract
sampcols=$(echo $(bcftools query -l $VCF | awk '{print 9+NR "\t" $1}' | grep -w -f $SAMPFILE | cut -f 1 ) | tr ' ' ,)
cols=$(echo "1-9,${sampcols}") # always need 1-9. samples start at 10

# Get header
bcftools view -h $VCF | grep "^##" > ${OUTPREFIX}.vcf
bcftools view -h $VCF | grep "^#CHROM" | cut -f${cols} >> ${OUTPREFIX}.vcf

# Get rest of data
gcloud storage cp $VCF - --billing-project=$GOOGLEPROJECT| bgzip --decompress | grep -v "^#" | cut -f${cols} >> ${OUTPREFIX}.vcf
bgzip ${OUTPREFIX}.vcf > ${OUTPREFIX}.vcf.gz && tabix -p vcf ${OUTPREFIX}.vcf.gz
exit 0

