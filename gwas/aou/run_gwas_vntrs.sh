#!/bin/bash

ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"



# Disease
samples="samples/EUR_WHITE.csv"
for chr in "chr16"; do
    aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
    for phenotype in "Malignant_neoplasm_of_the_skin"; do
    #for phenotype in "Hyperlipidemia"; do
    #for phenotype in "glaucoma"; do
            echo "Running gwas for $phenotype"
            python plot_genotype_phenotype.py --phenotype $phenotype \
                   --tr-vcf $aou_220k_imputed \
                   --samples $samples \
                   --outdir outputs/${chr} \
                   --annotations annotation_points.csv \
                   --binary
           : '
             vntr_plink_gwas_file="../../../rebased_cast-workflows/cast-workflows/tr-gwas/results_v3_eur_significant_in_v2/sara_${phenotype}_EUR_WHITE_gwas.tab"
            processed_vntr_gwas_file="data/plink/processed_gwas_${phenotype}.tab"
            cat $vntr_plink_gwas_file | grep "^#CHROM" | head -n 1 > $processed_vntr_gwas_file
            cat $vntr_plink_gwas_file | grep "ADD" >> $processed_vntr_gwas_file
            phecode="EM_239.1"
            snp_gwas_file="data/all_by_all/df_dump_${chr}_${phecode}.csv"
                ./aou_gwas.py --phenotype $phenotype \
                   --num-pcs 10 \
                   --method associaTR \
                   --tr-vcf $aou_220k_imputed \
                   --samples $samples \
                   --outdir outputs/${chr} \
                   --annotations annotation_points.csv \
                   --snp-gwas-file $snp_gwas_file \
                   --vntr-gwas-file $processed_vntr_gwas_file \
                   --plot-genotype-phenotype \
                   --binary \
                   --plot'
    done
done

exit 0

# Blood
samples="samples/passing_samples_v7.1.csv"
for chr in "chr1"; do
    #for phenotype in "red_blood_cell_distribution_width" "platelet_count" "mean_platelet_volume" "mean_corpuscular_hemoglobin" "cholesterol" "alkaline_phosphatase"; do
    for phenotype in "alkaline_phosphatase"; do
        #phecode="3035995"
        aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
        #snp_gwas_file="data/all_by_all/df_dump_${chr}_${phecode}.csv"
        echo "Running gwas for $phenotype"
        python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --norm quantile

       : " ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --plot-genotype-phenotype \
           --norm quantile \
           --plot
           #--snp-gwas-file $snp_gwas_file "
    done
done

exit 0


#most_significant_hit=$(tail -n +4 outputs/${chr}/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)
