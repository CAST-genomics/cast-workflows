#!/bin/bash
PHEN=$1
LOCI=$2
PFILE_PREFIX=$3

while IFS=$' ' read -r -a LOCI
do
    chr=${locus[0]}
    from=${locus[1]}
    to=${locus[2]}
    gsutil -q stat ${WORKSPACE_BUCKET}/nichole/plink/finemapping/${chr}_${from}_${to}_${phen}_str_plink.bed
    status=$?
    if [[ $status == 0 ]]; then
        echo "plink file exists"
    else
        plink2 --pfile "$PFILE_PREFIX" \
           --chr "$chr" \
           --from-bp "$from" \
           --to-bp "$to" \
           --make-pgen \
           --out "${chr}_${from}_${to}_${PHEN}_str_plink"     
    fi
done < $LOCI