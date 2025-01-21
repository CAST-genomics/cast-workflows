chr=$1
from=$2
to=$3
snppos=$4
gnomixchr=$5
gnomixstart=$6
gnomixend=$7
phenname=$8
phen=$9


gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phenname"
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.*" ./

gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phenname"
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.*" ./


samples=passing_samples_v7.1
snpvid=$(awk -v snppos="$snppos" '$4 == snppos {print $2}' "${chr}_${from}_${to}_${phenname}_plink.bim")

plink2 --bfile ${chr}_${from}_${to}_${phenname}_plink --snp $snpvid --make-bed --out $snppos
numsamples=$(wc -l < "${snppos}.fam")
python ../convert_bed_to_numpy.py ${snppos}.bed $numsamples 1
#above will write out ${snppos}.npy

numsamples=$(wc -l < "${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.fam")
numvars=$(wc -l < "${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.bim")
python ../convert_bed_to_numpy.py ${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.bed $numsamples $numvars

python get_snp_phenocovar.py "${samples}.csv" "$phen" "AOU_10_PCS.tsv" $snppos


python run_chr_variant_interaction_gwas.py "${samples}_${snppos}" ${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.npy ${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.bim ${gnomixchr}_${gnomixstart}_${gnomixend}_${phenname}_plink.fam 


#rm ${samples}_${snppos}
#rm $snppos.*
#rm ${chr}_${from}_${to}_${phenname}_plink.bed
#rm ${chr}_${from}_${to}_${phenname}_plink.bim
#rm ${chr}_${from}_${to}_${phenname}_plink.fam
