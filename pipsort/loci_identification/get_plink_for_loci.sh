phen=$1
loci=$2

while IFS=$' ' read -r -a locus
do
    chr=${locus[0]}
    from=${locus[1]}
    to=${locus[2]}
    gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phen}_plink.bed
    status=$?
    if [[ $status == 0 ]]; then
        echo "plink file exists"
    else
        python subset_hail_to_plink.py $chr $from $to $phen	    
    fi
done < $loci

