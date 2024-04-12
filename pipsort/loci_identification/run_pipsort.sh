chr=$1
from=$2
to=$3
sp=$4
fp=${sp#*.}
s1_samples=$5
s2_samples=$6
s1_gwas=$7
s2_gwas=$8
phen=$9 #e.g. ldl_cholesterol
#cojo_numsignals=$4
#cr_numsignals=$4

fsl_s1=0.0001
fsl_s2=0.0001
maf=0.01
min_num_snps=25
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/plink_bed/acaf_threshold.chr${chr}
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr${chr}

logfile=pipsort.log
results=results_ldl_${fp}.txt

#0. copy over samples files if needed
if [ -e "$s1_samples" ]; then
    echo "sample s1 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s1_samples}" ./
fi
if [ -e "$s2_samples" ]; then
    echo "sample file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s2_samples}" ./
fi

#1. copy over plink files if needed
plink_file=acaf_threshold.chr$chr.bed
if [ -e "$plink_file" ]; then
    echo "no need to copy plink data"
else
    echo ${plink_file_prefix}.*
    gsutil -u $GOOGLE_PROJECT -m cp ${plink_file_prefix}.* ./
fi
#!gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr21.fam ./


#2. get phenotype file for gwas
if [ -e "${phen}_phenocovar.csv" ]; then
    echo "no need to copy phen data"
else
    gsutil cp "${WORKSPACE_BUCKET}/phenotypes/${phen}_phenocovar.csv" ./
fi

exit 0

#3. subset to region I need
plink2 --bfile acaf_threshold.chr${chr} --maf $maf --chr $chr --from-bp $from --to-bp $to --keep $s1_samples --make-bed --out s1_data --pheno ${phen}_phenocovar.csv --prune
plink2 --bfile acaf_threshold.chr${chr} --maf $maf --chr $chr --from-bp $from --to-bp $to --keep $s2_samples --make-bed --out s2_data ${phen}_phenocovar.csv --prune




