chr=$1
from=$2
to=$3
snppos=$4
phenname=$5
phen=$6

gnomix="${WORKSPACE_BUCKET}/gnomix/outputs/gnomix-chr${chr}.msp"

# Check if the file exists in the bucket
if [ -f "gnomix-chr${chr}.msp" ]; then
    echo "File ${gnomix} exists."
else
    gsutil cp "${gnomix}" .
fi

awk -v pos="$snppos" '$2 <= pos && pos <= $3' gnomix-chr${chr}.msp > temp
sed -n '2p' gnomix-chr${chr}.msp > temp2

cat temp2 > region.txt
cat temp >> region.txt

rm temp
rm temp2

python gnomix_to_local_ancestry.py region.txt region_lancestry.tsv eur

gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phenname"
    exit 1
fi
gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.*" ./


lancestry_codes=(0 1 2)
samples_list=("AFR_BLACK" "EUR_WHITE" "NOT_AFR_BLACK" "passing_samples_v7.1")
snpvid=$(awk -v snppos="$snppos" '$4 == snppos {print $2}' "${chr}_${from}_${to}_${phenname}_plink.bim")

resultsfile=${chr}_${from}_${to}_${snpvid}_results.csv
echo "group,interval,n,beta,se,p" > $resultsfile

for lancestry_code in "${lancestry_codes[@]}"; do
  for samples in "${samples_list[@]}"; do
    echo "Processing lancestry_code=${lancestry_code}, samples=${samples}..."

    # Run the Python script
    python get_local_ancestry_phencovar.py "${samples}.csv" "$phen" "AOU_10_PCS.tsv" "region_lancestry.tsv" "${lancestry_code}" "eur"

    # Get SNP ID

    # Run plink2
    plink2 --bfile "${chr}_${from}_${to}_${phenname}_plink" \
           --snp "$snpvid" \
           --linear hide-covar \
           --pheno "${samples}_eur_${lancestry_code}" \
           --pheno-name phenotype \
           --covar "${samples}_eur_${lancestry_code}" \
           --covar-name age-pc10 \
           --covar-variance-standardize \
           --vif 2000

    label="${lancestry_code},$samples"
    # Extract fields from the output file
    if [[ -f plink2.phenotype.glm.linear ]]; then
      awk -v label="$label" 'NR > 1 {print label "," $11 "," $12 "," $13 "," $15}' plink2.phenotype.glm.linear >> $resultsfile
      rm plink2.phenotype.glm.linear
    else
      echo "plink2.phenotype.glm.linear not found for lancestry_code=${lancestry_code}, samples=${samples}"
    fi
    rm ${samples}_eur_${lancestry_code}
  done
done

rm region_lancestry.tsv
rm ${chr}_${from}_${to}_${phenname}_plink.bed
rm ${chr}_${from}_${to}_${phenname}_plink.bim
rm ${chr}_${from}_${to}_${phenname}_plink.fam
