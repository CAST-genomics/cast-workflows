#!/bin/bash

echo "" > summary_gwas.txt
#aou_220k_imputed="/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/vntr_workspace/sara_main/cast-workflows/imputation/wdl/batch_data/merged_samples.sorted.vcf.gz"
aou_220k_imputed="merged_samples.sorted_with_vid_ru.vcf.gz"
lrwgs_data="../../../../imputation/cast-workflows/imputation/wdl/vntr_reference/merged_samples.sorted_clean.vcf"


#samples="samples/passing_samples_v7.1.csv"
#samples="samples/EUR_WHITE.csv"
#samples="samples/AFR_BLACK.csv"


# Running associatr
echo "running gwas for imputed calls"
samples="samples/EUR_WHITE.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --annotations annotations_acan.txt \
	      --samples $samples \
	      --is-imputed \
	      --norm-by-sex \
	      --norm quantile \
	      --plot
exit 0
samples="samples/AFR_BLACK.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --annotations annotations_acan.txt \
	      --samples $samples \
	      --is-imputed \
	      --norm-by-sex \
	      --norm quantile \
	      --plot
exit 0
echo "running gwas for advntr calls"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $lrwgs_data \
	      --annotations annotations_acan.txt \
	      --norm quantile \
	      --samples $samples \
	      --norm-by-sex \
	      --annotations annotations_acan.txt \
	      --plot
#	      --tr-vcf merged_samples.vcf \
# Running Hail
#./aou_gwas.py --phenotype height \
#	      --num-pcs 10 \
#	      --region chr15:77855424-99857434 \
#	      --method hail \
#	      --norm quantile \
#	      --norm-by-sex \
#	      --plot
	      #--region chr15:83855424-93857434 \ #5Mbp window
	      #--region chr15:87855424-89857434 \ # 1Mbp window
	      #--region chr15:88000000-90000000 \
	      #--region chr15:88855400-88857500 \
	      #--region chr15:88855424-88857434
	      #--tr-vcf ${WORKSPACE_BUCKET}/saraj/batch_genotyping/run_8/merged_outputs/run_merge_advntr/d06c3f27-b7b7-41ad-9cfb-8e7161c43fce/call-merge_outputs/merged_samples.vcf

exit 0

for phenotype in $(tail -n +2 phenotypes_manifest.csv | cut -d, -f1); do
	./aou_gwas.py --phenotype $phenotype \
              --num-pcs 10 \
              --method associaTR \
              --tr-vcf merged_samples.vcf \
              --norm zscore \
              --plot
      most_significant_hit=$(tail -n +4 outputs/${phenotype}_associaTR_ALL.gwas.tab | cut -f6 | sort  | awk NF | head -n 1)
      echo "most_significant_hit for phenotype $phenotype is $most_significant_hit" >> summary_gwas.txt
      echo "----------- most_significant_hit for phenotype $phenotype is $most_significant_hit"
done
