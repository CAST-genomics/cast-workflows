chr=$1
snpposhg19=$2

gnomix="${WORKSPACE_BUCKET}/gnomix/outputs/gnomix-chr${chr}.msp"

# Check if the file exists in the bucket
if [ -f "gnomix-chr${chr}.msp" ]; then
    echo "File ${gnomix} exists."
else
    gsutil cp "${gnomix}" .
fi

awk -v pos="$snpposhg19" '$2 <= pos && pos <= $3' gnomix-chr${chr}.msp > temp
sed -n '2p' gnomix-chr${chr}.msp > temp2

cat temp2 > region.txt
cat temp >> region.txt

rm temp
rm temp2

python gnomix_to_local_ancestry.py region.txt region_lancestry.tsv nat 2
