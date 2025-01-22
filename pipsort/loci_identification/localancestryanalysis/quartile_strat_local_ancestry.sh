chr=$1
from=$2
to=$3
snppos=$4
phenname=$5
phen=$6


gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phenname"
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.*" ./


samples_list=("EUR_Q1_samples" "EUR_Q2_samples" "EUR_Q3_samples" "EUR_Q4_samples")
snpvid=$(awk -v snppos="$snppos" '$4 == snppos {print $2}' "${chr}_${from}_${to}_${phenname}_plink.bim")

resultsfile=${chr}_${from}_${to}_${snpvid}_${phenname}_results.csv
echo "group,n,beta,se,p" > $resultsfile

for samples in "${samples_list[@]}"; do

  # Run the Python script
  python get_phenocovar.py "${samples}.csv" "$phen" "AOU_10_PCS.tsv" ${samples}_${snppos}.csv

  # Run plink2
  plink2 --bfile "${chr}_${from}_${to}_${phenname}_plink" \
         --snp "$snpvid" \
         --linear hide-covar \
         --pheno "${samples}_${snppos}.csv" \
         --pheno-name phenotype \
         --covar "${samples}_${snppos}.csv" \
         --covar-name age-pc10 \
         --covar-variance-standardize \
         --vif 2000

  label="$samples"
  # Extract fields from the output file
  if [[ -f plink2.phenotype.glm.linear ]]; then
    awk -v label="$label" 'NR > 1 {print label "," $11 "," $12 "," $13 "," $15}' plink2.phenotype.glm.linear >> $resultsfile
    rm plink2.phenotype.glm.linear
  else
    echo "plink2.phenotype.glm.linear not found for samples=${samples}"
  fi
  rm ${samples}_${snppos}.csv
done

rm ${chr}_${from}_${to}_${phenname}_plink.bed
rm ${chr}_${from}_${to}_${phenname}_plink.bim
rm ${chr}_${from}_${to}_${phenname}_plink.fam
