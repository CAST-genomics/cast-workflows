#!/bin/bash

ref="data/lrwgs_p_g_polymorphic_vntrs_sr_6_ml_95.sorted.vcf.gz"

# Disease

samples="samples/EUR_WHITE.csv"
samples_prefix="EUR"

chr="chr1"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Abnormal_intraocular_pressure"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done

exit 0

chr="chr16"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Actinic_keratosis" "Malignant_neoplasm_of_the_skin" "Skin_changes_due_to_chronic_exposure_to_nonionizing_radiation" "Skin_changes_due_to_solar_radiation" "Melanomas_of_skin" "Keratinocyte_carcinoma" "Basal_cell_carcinoma"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done

chr="chr2"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Hypercholesterolemia" "Hyperlipidemia"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


exit 0

chr="chr12"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Open_angle_glaucoma"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done

exit 0

           
             vntr_plink_gwas_file="../../../rebased_cast-workflows/cast-workflows/tr-gwas/results_v3_eur_significant_in_v2/sara_${phenotype}_EUR_WHITE_gwas.tab"
            processed_vntr_gwas_file="data/plink/processed_gwas_${phenotype}.tab"
            cat $vntr_plink_gwas_file | grep "^#CHROM" | head -n 1 > $processed_vntr_gwas_file
            cat $vntr_plink_gwas_file | grep "ADD" >> $processed_vntr_gwas_file
            snp_gwas_file="data/all_by_all/df_dump_${chr}_${phenotype}.csv"
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
                   --plot
    


# Blood
samples="samples/passing_samples_v7.1.csv"
for chr in "chr2"; do
    #for phenotype in "red_blood_cell_distribution_width" "platelet_count" "mean_platelet_volume" "mean_corpuscular_hemoglobin" "cholesterol" "alkaline_phosphatase"; do
    for phenotype in "cholesterol"; do
        phecode=$phenotype
        aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
        snp_gwas_file="data/all_by_all/df_dump_${chr}_${phecode}.csv"
        echo "Running gwas for $phenotype"
        : "python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --norm zscore"

        ./aou_gwas.py --phenotype $phenotype \
           --num-pcs 10 \
           --method associaTR \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --annotations annotation_points.csv \
           --plot-genotype-phenotype \
           --snp-gwas-file $snp_gwas_file \
           --norm zscore \
           --plot
    done
done



exit 0

exit 0


#most_significant_hit=$(tail -n +4 outputs/${chr}/${phenotype}_associaTR_${samples_prefix}.gwas.tab | cut -f6 | sort -g  | awk NF | head -n 1)


exit 0
chr="chr6"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Actinic_keratosis" "Skin_changes_due_to_chronic_exposure_to_nonionizing_radiation"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done
chr="chr7"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Basal_cell_carcinoma" "Actinic_keratosis"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


chr="chr20"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Actinic_keratosis" "Skin_changes_due_to_chronic_exposure_to_nonionizing_radiation"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


chr="chr21"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Skin_changes_due_to_solar_radiation"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


chr="chr18"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Chronic_thyroiditis" "Hashimoto_thyroiditis__Chronic_lymphocytic_thyroiditis_" "Rheumatic_fever_and_chronic_rheumatic_heart_diseases"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


samples="samples/AMR_HISPANIC.csv"
samples_prefix="AMR"

chr="chr6"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Skin_changes_due_to_chronic_exposure_to_nonionizing_radiation"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


samples="samples/AFR_BLACK.csv"
samples_prefix="AFR"



chr="chr10"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Melanomas_of_skin"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done

exit 0

chr="chr11"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Basal_cell_carcinoma"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done


chr="chr14"
aou_220k_imputed="../../imputation/wdl/data/imputed_${chr}.annotated.rh.vcf.gz"
for phenotype in "Basal_cell_carcinoma"; do
    output="out_${chr}_${phenotype}_${samples_prefix}.txt"
    rm -f $output
    echo "Running gwas for $phenotype"
    python plot_genotype_phenotype.py --phenotype $phenotype \
           --tr-vcf $aou_220k_imputed \
           --samples $samples \
           --outdir outputs/${chr} \
           --binary \
           --annotations annotation_points.csv >> $output
done
