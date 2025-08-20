#full pipeline for running plink on a region
chr=$1
from=$2
to=$3
phenname=$4
samples=$5 #options are passing_samples_v7.1.csv, EUR_WHITE.csv, AFR_BLACK.csv, hispanic
norm=False
if [ $# -ge 6 ]; then
  norm=$6
fi

scripts=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/

# get plink files
gsutil -q stat ${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.bed
status=$?
if [[ $status == 0 ]]; then
    echo "plink file exists"
else
    echo "plink file does not exist $chr $to $from $phenname"
    exit 1
fi
#gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_${phenname}_plink.*" ./


# get pcs
if [ ! -f "AOU_10_PCS.tsv" ]; then
    gsutil cp "${WORKSPACE_BUCKET}/tmirmira/data/covariates/v7/AOU_10_PCS.tsv" .
fi

# get all samples
if [ ! -f "${samples}" ]; then
    gsutil cp "${WORKSPACE_BUCKET}/samples/${samples}" .
fi


# get phenocovar. For now, I will give it a file TODO
phenocovar=${phenname}_phenocovar.csv


python $scripts/localancestryanalysis/get_phenocovar.py $samples ${phenname}_phenocovar.csv "AOU_10_PCS.tsv" ${phenname}_plink.csv $norm


# Run plink2
plink2 --bfile "${chr}_${from}_${to}_${phenname}_plink" \
       --linear hide-covar \
       --pheno "${phenname}_plink.csv" \
       --pheno-name phenotype \
       --covar "${phenname}_plink.csv" \
       --covar-name age-pc10 \
       --covar-variance-standardize \
       --maf 0.01 \
       --vif 2000

basename="${samples%.csv}"
python plot_gwas_plink.py plink2.phenotype.glm.linear ${basename}_norm${norm}.png


