chr=$1
from=$2
to=$3
snppos=$4
snpposhg19=$5
phenname=$6
phen=$7


gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phenname"
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.*" ./

samples=passing_samples_v7.1
snpvid=$(awk -v snppos="$snppos" '$4 == snppos {print $2}' "${chr}_${from}_${to}_${phenname}_plink.bim")

plink2 --bfile ${chr}_${from}_${to}_${phenname}_plink --snp $snpvid --make-bed --out $snppos
numsamples=$(wc -l < "${snppos}.fam")
python ../convert_bed_to_numpy.py ${snppos}.bed $numsamples 1
#above will write out ${snppos}.npy

python get_snp_phenocovar.py "${samples}.csv" "$phen" "AOU_10_PCS.tsv" $snppos


local_file="chr${chr}_local_ancestry_0.csv"
if [[ ! -f "$local_file" ]]; then
	echo "File $local_file not found locally. Copying from bucket..."
	gsutil cp "${WORKSPACE_BUCKET}/pipsort/localancestry/eur/${local_file}" .
fi
output_file="chr${chr}_${snppos}_interaction_results_no_het.csv"
python run_interaction_no_het.py "${samples}_${snppos}" "$local_file" $snpposhg19

rm ${samples}_${snppos}
rm $snppos.*
rm ${chr}_${from}_${to}_${phenname}_plink.bed
rm ${chr}_${from}_${to}_${phenname}_plink.bim
rm ${chr}_${from}_${to}_${phenname}_plink.fam
