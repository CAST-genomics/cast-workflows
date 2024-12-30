chr=$1
from=$2
to=$3
snppos=$4
phenname=$5
phen=$6

logfile="alllogsgwasinteraction"

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

samples=passing_samples_v7.1
snpvid=$(awk -v snppos="$snppos" '$4 == snppos {print $2}' "${chr}_${from}_${to}_${phenname}_plink.bim")

plink2 --bfile ${chr}_${from}_${to}_${phenname}_plink --snp $snpvid --make-bed --out $snppos
numsamples=$(wc -l < "${snppos}.fam")
python ../convert_bed_to_numpy.py ${snppos}.bed $numsamples 1
#above will write out ${snppos}.npy

python get_snp_phenocovar.py "${samples}.csv" "$phen" "AOU_10_PCS.tsv" "region_lancestry.tsv" $snppos

final_results="${chr}_${snppos}_interaction_results.csv"
> "$final_results" # Clear or create the file

for chriter in {21..22}; do
	echo "$chriter INTERACTION GWAS"
	local_file="chr${chriter}_local_ancestry_0.csv"
	if [[ ! -f "$local_file" ]]; then
		echo "File $local_file not found locally. Copying from bucket..."
		gsutil cp "${WORKSPACE_BUCKET}/pipsort/localancestry/eur/${local_file}" .
	fi
	output_file="chr${chriter}_${snppos}_interaction_results.csv"
	python run_chr_interaction_gwas.py "${samples}_${snppos}" "$local_file"

	if [[ -f "$output_file" ]]; then
		echo "Appending results from $output_file to $final_results..."
		cat "$output_file" >> "$final_results"
		rm $output_file
	else
		echo "Warning: Output file $output_file not found." >> $logfile
	fi
done

rm ${samples}_${snppos}
rm $snppos.*
rm region.txt
rm region_lancestry.tsv
rm ${chr}_${from}_${to}_${phenname}_plink.bed
rm ${chr}_${from}_${to}_${phenname}_plink.bim
rm ${chr}_${from}_${to}_${phenname}_plink.fam
