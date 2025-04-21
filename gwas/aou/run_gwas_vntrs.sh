#!/bin/bash

ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"

# Disease
samples="samples/EUR_WHITE.csv"
for chr in "chr2"; do # "chr6" "chr8" "chr9" "chr22"; do
    aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
    for phenotype in "Hyperlipidemia"; do
            echo "Running gwas for $phenotype"
            ./aou_gwas.py --phenotype $phenotype \
                   --num-pcs 10 \
                   --method associaTR \
                   --tr-vcf $aou_220k_imputed \
                   --samples $samples \
                   --outdir outputs/${chr} \
                   --annotations annotation_points.csv \
                   --plot-genotype-phenotype \
                   --binary \
                   --plot
    done
done

exit 0
# Blood
#samples="samples/passing_samples_v7.1.csv"
#samples="samples/AFR_BLACK.csv"
#samples="samples/EUR_WHITE.csv"
for chr in "chr2"; do # "chr11" "chr16"
    #for phenotype in "red_blood_cell_distribution_width" "platelet_count" "mean_platelet_volume" "mean_corpuscular_hemoglobin"; do
    for phenotype in "cholesterol"; do
        aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
        snp_gwas_file="data/all_by_all/df_dump_${chr}_${phenotype}.csv"
        echo "Running gwas for $phenotype"
        ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --plot-genotype-phenotype \
           --norm quantile \
           --plot
           #--snp-gwas-file $snp_gwas_file \
    done
done





#most_significant_hit=$(tail -n +4 outputs/${chr}/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)
#echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> $summary
#echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
