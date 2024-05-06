chr=$1
from=$2
to=$3
num_sig_snps_in_loci=$4
s1_samples_file=$5
s1_samples="${s1_samples_file%.*}"
s2_samples_file=$6
s2_samples="${s2_samples_file%.*}"
phen=$7 #e.g. ldl_cholesterol

max_num_cr=15
fsl_s1=0.0001
fsl_s2=0.0001
maf=0.01
min_num_snps=25
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/plink_bed/acaf_threshold.chr${chr}
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr${chr}

plink_loc=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification
scripts=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification
common=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/pipsort_inputs_ldl_cholesterol


cd /home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/pipsort_inputs_${phen}
logfile=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/pipsort_inputs_${phen}/pipsort.log
cr_info=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/pipsort_inputs_${phen}/loci_cr_info


if [ -d "${chr}_${from}_${to}" ]; then
        echo "delete ${chr}_${from}_${to}"
        rm -r "${chr}_${from}_${to}"
fi

mkdir ${chr}_${from}_${to}
cd ${chr}_${from}_${to}


#0. copy over samples files if needed
if [ -e "$s1_samples_file" ]; then
    echo "sample s1 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s1_samples_file}" ./
fi
if [ -e "$s2_samples_file" ]; then
    echo "sample s2 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s2_samples_file}" ./
fi

#1. get plink files
gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phen}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phen" >> $logfile
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phen}_plink.*" ./


#2. get phenotype file for gwas
if [ -e "${phen}_phenocovar.csv" ]; then
    echo "no need to copy phen data"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/phenotypes/${phen}_phenocovar.csv" ./
fi
python $scripts/convert_phen_to_plink_format.py ${phen}_phenocovar.csv ${phen}_plink_format.tab

#3. get gwas data and add rsid to gwas file
gwas_file_s1_pre=${phen}_hail_${s1_samples}
gwas_file_s2_pre=${phen}_hail_${s2_samples}

if [ -e "$common/${gwas_file_s1_pre}.gwas.tab" ]; then
    echo "no need to copy gwas1"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/gwas/${gwas_file_s1_pre}.gwas.tab" ./
fi
if [ -e "$common/${gwas_file_s2_pre}.gwas.tab" ]; then
    echo "no need to copy gwas2"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/gwas/${gwas_file_s2_pre}.gwas.tab" ./
fi


python $scripts/subset_gwas_to_loci.py $chr $from $to $common/${gwas_file_s1_pre}.gwas.tab ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab
python $scripts/subset_gwas_to_loci.py $chr $from $to $common/${gwas_file_s2_pre}.gwas.tab ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab
python $scripts/add_rsid_col.py ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab
python $scripts/add_rsid_col.py ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab


#3. subset to snps and samples I need
python $scripts/convert_samples_to_plink_format.py $s1_samples_file plink_${s1_samples_file}
python $scripts/convert_samples_to_plink_format.py $s2_samples_file plink_${s2_samples_file}
python $scripts/extract_all_snps.py ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab union_snps.txt

$plink_loc/plink2 --bfile ${chr}_${from}_${to}_${phen}_plink --chr $chr --keep plink_$s1_samples_file --make-bed --out s1_data --pheno ${phen}_plink_format.tab --prune --extract union_snps.txt
$plink_loc/plink2 --bfile ${chr}_${from}_${to}_${phen}_plink --chr $chr --keep plink_$s2_samples_file --make-bed --out s2_data --pheno ${phen}_plink_format.tab --prune --extract union_snps.txt

python $scripts/sort_gwas_results.py --infile ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab --pos_col pos --rsid_col rsid
python $scripts/sort_gwas_results.py --infile ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab --pos_col pos --rsid_col rsid

#cond reg analysis
bash $scripts/single_study_cond_reg_with_input_data.sh s1_data ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab ${phen}_phenocovar.csv 0.0001 s1_lead_snps $max_num_cr
bash $scripts/single_study_cond_reg_with_input_data.sh s2_data ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab ${phen}_phenocovar.csv 0.0001 s2_lead_snps $max_num_cr
if [ -e "s1_lead_snps" ]; then
    echo "s1_lead_snps exists"
else
    # Create the file if it doesn't exist
    touch "s1_lead_snps"
fi
if [ -e "s2_lead_snps" ]; then
    echo "s2_lead_snps exists"
else
    # Create the file if it doesn't exist
    touch "s2_lead_snps"
fi
a=($(wc -l s1_lead_snps))
num_lead_snps_s1=${a[0]}
a=($(wc -l s2_lead_snps))
num_lead_snps_s2=${a[0]}


#exclude high p vals, outputs s1_snps.txt and s2_snps.txt and s1_gwas_temp.txt and s2_gwas_temp.txt
python $scripts/remove_high_pvals_not_common.py --gwas1 ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab --gwas2 ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab --fsl1 $fsl_s1 --fsl2 $fsl_s2 --pval_col p_value --rsid_col rsid

mv ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab s1_orig_gwas
mv ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab s2_orig_gwas

python $scripts/get_gwas_in_mscaviar_style.py --infile s1_gwas_temp.txt --outfile s1_gwas.txt --chr_col chrom --pos_col pos --rsid_col rsid --ref_col n --alt_col locus --beta_col beta --se_col standard_error
python $scripts/get_gwas_in_mscaviar_style.py --infile s2_gwas_temp.txt --outfile s2_gwas.txt --chr_col chrom --pos_col pos --rsid_col rsid --ref_col n --alt_col locus --beta_col beta --se_col standard_error
#rm s1_gwas_temp.txt
#rm s2_gwas_temp.txt

#LD
$plink_loc/plink2 --bfile s1_data  --extract s1_snps.txt --make-bed --out s1_remaining_snps
$plink_loc/plink2 --bfile s2_data  --extract s2_snps.txt --make-bed --out s2_remaining_snps
$plink_loc/plink2 --bfile s1_remaining_snps --r-unphased square
mv plink2.unphased.vcor1 s1.ld
$plink_loc/plink2 --bfile s2_remaining_snps --r-unphased square
mv plink2.unphased.vcor1 s2.ld

#AFREQ
$plink_loc/plink2 --bfile s1_remaining_snps --freq
mv plink2.afreq s1.afreq
$plink_loc/plink2 --bfile s2_remaining_snps --freq
mv plink2.afreq s2.afreq


mscav=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/MsCAVIAR

utils=$mscav/utils
python $utils/format_sumstats.py --infile s1_gwas.txt --outfile f_s1.txt --chromosome CHROMOSOME --bp BP --snp_id SNP_ID --ref_allele REF_ALLELE --alt_allele ALT_ALLELE --beta BETA --se SE
python $utils/format_sumstats.py --infile s2_gwas.txt --outfile f_s2.txt --chromosome CHROMOSOME --bp BP --snp_id SNP_ID --ref_allele REF_ALLELE --alt_allele ALT_ALLELE --beta BETA --se SE

python $scripts/format_zscores.py f_s1.txt s1_final.zscores
python $scripts/format_zscores.py f_s2.txt s2_final.zscores


a=($(wc -l s1_final.zscores))
num_snps_s1=${a[0]}
if test $num_snps_s1 -lt $min_num_snps
then
    echo "Study 1 has ${num_snps_s1} SNPs"
    echo "*** $1 $2 $3 FAILED***: study 1 too few snps" >> $logfile
    cd ..
    rm -r "${chr}_${from}_${to}"
    #rm s1*
    #rm s2*
    #rm $s1_samples_file
    #rm $s2_samples_file
    #rm ${phen}_phenocovar.csv
    #rm plink*
    #rm f_s1.txt
    #rm f_s2.txt
    exit 0
fi

a=($(wc -l s2_final.zscores))
num_snps_s2=${a[0]}
if test $num_snps_s2 -lt $min_num_snps
then
    echo "Study 2 has ${num_snps_s2} SNPs"
    echo "*** $1 $2 $3 FAILED***: study 2 too few snps" >> $logfile
    cd ..
    rm -r "${chr}_${from}_${to}"
    #rm s1*
    #rm s2*
    #rm $s1_samples_file
    #rm $s2_samples_file
    #rm ${phen}_phenocovar.csv
    #rm plink*
    #rm f_s1.txt
    #rm f_s2.txt
    exit 0
fi

echo "$chr $from $to $num_sig_snps_in_loci $num_lead_snps_s1 $num_lead_snps_s2" >> $cr_info

echo $'f_s1.txt\nf_s2.txt' > processedfiles.txt
python $scripts/extract_all_snps_rsid.py processedfiles.txt snp_map



#cleanup all files
rm s1_lead_snps
rm s2_lead_snps
rm s1_snps.txt
rm s2_snps.txt
rm union_snps.txt
rm f_s1.txt
rm f_s2.txt
rm *.tab
rm $s1_samples_file
rm $s2_samples_file
rm *.bed
rm *.fam
rm *.bim
#rm snp_map
rm processedfiles.txt
rm plink*
rm *.log
rm ${phen}_phenocovar.csv
