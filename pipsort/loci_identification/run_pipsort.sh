chr=$1
from=$2
to=$3
sp=$4
fp=${sp#*.}
s1_samples=$5
s2_samples=$6
s1_gwas=$7
s2_gwas=$8
#cojo_numsignals=$4
#cr_numsignals=$4

fsl_s1=0.0001
fsl_s2=0.0001
maf=0.01
min_num_snps=25
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/plink_bed/acaf_threshold.chr${chr}

logfile=pipsort.log
results=results_ldl_${fp}.txt

#1. copy over plink files if needed
plink_file=acaf_threshold.chr$chr.bed
if [ -e "$plink_file" ]; then
    echo "no need to copy plink data"
else
    gsutil -u $GOOGLE_PROJECT -m cp $plink_file_prefix.*
fi

#2. get phenotype file for gwas


#3. subset to region I need
plink2 --bfile acaf_threshold.chr${chr} --maf $maf --chr $chr --from-bp $from --to-bp $to --keep $s1_samples --make-bed --out s1_data
plink2 --bfile acaf_threshold.chr${chr} --maf $maf --chr $chr --from-bp $from --to-bp $to --keep $s2_samples --make-bed --out s2_data




