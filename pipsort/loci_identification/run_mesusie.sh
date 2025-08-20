chr=$1
from=$2
to=$3
s1_samples_file=$4
s1_samples="${s1_samples_file%.*}"
s2_samples_file=$5
s2_samples="${s2_samples_file%.*}"
phen=$6 #e.g. ldl_cholesterol

maf=0.01
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/plink_bed/acaf_threshold.chr${chr}
plink_file_prefix=gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr${chr}

plink_loc=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification
plink_loc=/usr/bin
scripts=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification
common=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/mesusie

if [ ! -d "$common" ]; then
        mkdir -p "$common"
        echo "Directory '$common' created."
    else
        echo "Directory '$common' already exists."
    fi

cd $common

logfile=$common/mesusie.log
if [ ! -f "$logfile" ]; then
        touch "$logfile"
        echo "File '$logfile' created."
    else
        echo "File '$logfile' already exists."
    fi


if [ -d "${chr}_${from}_${to}" ]; then
        echo "delete ${chr}_${from}_${to}"
        rm -r "${chr}_${from}_${to}"
fi

mkdir ${chr}_${from}_${to}
cd ${chr}_${from}_${to}


#0. copy over samples files if needed
if [ -e "$common/$s1_samples_file" ]; then
    echo "sample s1 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s1_samples_file}" $common
fi
if [ -e "$common/$s2_samples_file" ]; then
    echo "sample s2 file already exists"
else
    gsutil cp "${WORKSPACE_BUCKET}/samples/${s2_samples_file}" $common
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
if [ -e "$common/${phen}_phenocovar.csv" ]; then
    echo "no need to copy phen data"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/phenotypes/${phen}_phenocovar.csv" $common
fi
python $scripts/convert_phen_to_plink_format.py $common/${phen}_phenocovar.csv ${phen}_plink_format.tab

#3. get gwas data and add rsid to gwas file
gwas_file_s1_pre=${phen}_hail_${s1_samples}
gwas_file_s2_pre=${phen}_hail_${s2_samples}

if [ -e "$common/${gwas_file_s1_pre}.gwas.tab" ]; then
    echo "no need to copy gwas1"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/gwas/${gwas_file_s1_pre}.gwas.tab" $common
fi
if [ -e "$common/${gwas_file_s2_pre}.gwas.tab" ]; then
    echo "no need to copy gwas2"
else
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/gwas/${gwas_file_s2_pre}.gwas.tab" $common
fi


python $scripts/subset_gwas_to_loci.py $chr $from $to $common/${gwas_file_s1_pre}.gwas.tab ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab
python $scripts/subset_gwas_to_loci.py $chr $from $to $common/${gwas_file_s2_pre}.gwas.tab ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab
python $scripts/add_rsid_col.py ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab
python $scripts/add_rsid_col.py ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab


#3. subset to snps and samples I need
python $scripts/convert_samples_to_plink_format.py $common/$s1_samples_file plink_${s1_samples_file}
python $scripts/convert_samples_to_plink_format.py $common/$s2_samples_file plink_${s2_samples_file}

python $scripts/sort_gwas_results.py --infile ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab --pos_col pos --rsid_col rsid
python $scripts/sort_gwas_results.py --infile ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab --pos_col pos --rsid_col rsid

fsl=1
python $scripts/remove_high_pvals_common.py --gwas1 ${gwas_file_s1_pre}_${chr}_${from}_${to}.gwas.tab --gwas2 ${gwas_file_s2_pre}_${chr}_${from}_${to}.gwas.tab --fsl $fsl --pval_col p_value --rsid_col rsid
#creates s1_snps.txt, s2_snps.txt s1_gwas_temp.txt, s2_gwas_temp.txt

$plink_loc/plink2 --bfile ${chr}_${from}_${to}_${phen}_plink --chr $chr --keep plink_$s1_samples_file --make-bed --out s1_data --pheno ${phen}_plink_format.tab --prune --extract s1_snps.txt
$plink_loc/plink2 --bfile ${chr}_${from}_${to}_${phen}_plink --chr $chr --keep plink_$s2_samples_file --make-bed --out s2_data --pheno ${phen}_plink_format.tab --prune --extract s2_snps.txt

#LD
$plink_loc/plink2 --bfile s1_data --r-unphased square ref-based
mv plink2.unphased.vcor1 s1.ld
$plink_loc/plink2 --bfile s2_data --r-unphased square ref-based
mv plink2.unphased.vcor1 s2.ld

ancestry_weight="c(3,3,1)"
nc=3
Rscript $scripts/run_mesusie.R ./mesusieoutput $nc $ancestry_weight ./mesusiecsinfo

#cleanup all files
rm f_s1.txt
rm f_s2.txt
rm *.tab
rm $s1_samples_file
rm $s2_samples_file
rm *.bed
rm *.fam
#rm *.bim
#rm snp_map
rm processedfiles.txt
rm plink*
rm *.log
rm ${phen}_phenocovar.csv
