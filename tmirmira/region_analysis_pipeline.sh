chr=$1
from=$2
to=$3
phenname=$4

scripts=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/


# get all samples
if [ ! -f "${phenname}_phenocovar.csv" ]; then
    gsutil cp "${WORKSPACE_BUCKET}/phenotypes/${phenname}_phenocovar.csv" .
fi

# get plink files
echo ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $from $to $phenname"
    exit 1
fi
if [ ! -f "${chr}_${from}_${to}_${phenname}_plink_maf001.bed" ]; then
    gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.*" .
    cat ${phenname}_phenocovar.csv | cut -f1 > sample_ids_only 
    python convert_samples_to_plink_format.py sample_ids_only sample_ids_plink_${phenname}.txt
    plink2 --bfile ${chr}_${from}_${to}_${phenname}_plink \
	    --maf 0.001 \
	    --keep sample_ids_plink_${phenname}.txt \
	    --make-bed \
	    --out ${chr}_${from}_${to}_${phenname}_plink_maf001
fi
