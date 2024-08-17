chr=$1
from=$2
to=$3

mkdir -p ${chr}_${from}_${to}_cluster
cd ${chr}_${from}_${to}_cluster

scripts=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/


gsutil cp "${WORKSPACE_BUCKET}/pipsort/plink/${chr}_${from}_${to}_indep_regions_plink.*" ./

maf=0.01

plink2 --bfile ${chr}_${from}_${to}_indep_regions_plink --maf ${maf} --make-bed --out ${chr}_${from}_${to}_indep_regions_plink_maf

a=($(wc -l ${chr}_${from}_${to}_indep_regions_plink_maf.bim))
num_variants=${a[0]}

a=($(wc -l ${chr}_${from}_${to}_indep_regions_plink_maf.fam))
num_samples=${a[0]}

python $scripts/convert_bed_to_numpy.py ${chr}_${from}_${to}_indep_regions_plink_maf.bed $num_samples $num_variants

python $scripts/cluster_region.py ${chr}_${from}_${to}_indep_regions_plink_maf.npy ${chr}_${from}_${to}_indep_regions_plink_maf

gsutil cp ${chr}_${from}_${to}_indep_regions_plink_maf_kmeans_labels.npy "gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/pipsort/indep_regions/AFR"
gsutil cp ${chr}_${from}_${to}_indep_regions_plink_maf_pcs.npy "gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/pipsort/indep_regions/AFR"

touch ${chr}_${from}_${to}_cluster_counts

python $scripts/process_clustering_results.py ${chr}_${from}_${to}_indep_regions_plink_maf.fam ${chr}_${from}_${to}_indep_regions_plink_maf_kmeans_labels.npy ${chr}_${from}_${to}_indep_regions_plink_maf_pcs.npy ${chr}_${from}_${to}_cluster_counts

