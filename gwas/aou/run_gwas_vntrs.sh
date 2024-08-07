#!/bin/bash

echo "" > summary_gwas.txt
#aou_100_imputed="/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/vntr_workspace/nichole_imputation/cast-workflows/imputation/wdl/data/output_chr15_10mb_aou_275_lrwgs.vcf"
#aou_10k_imputed="../../../../nichole_imputation/cast-workflows/imputation/wdl/data/output_chr15_acan_only_aou_10k_samples_imputed_rh.vcf"
aou_10k_imputed="../../imputation/wdl/data/output_chr15_acan_50mbp_aou_10k_srwgs_samples_output.vcf.gz"
aou_10k_batch_merged="../../imputation/wdl/batch_data/vntr_from_merged_samples_sorted.vcf.gz"
aou_50k_imputed="../../imputation/wdl/data/output_chr15_acan_50mbp_aou_50k_srwgs_samples_output.sorted.vcf.gz"
aou_1000_imputed_batch="/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/vntr_workspace/sara_main/cast-workflows/imputation/wdl/batch_data/batch_set1_output.sorted.vcf.gz"
#aou_220k_imputed="/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/vntr_workspace/sara_main/cast-workflows/imputation/wdl/batch_data/merged_samples.sorted.vcf.gz"
aou_220k_imputed="merged_samples.sorted_with_vid_ru.vcf.gz"
lrwgs_data="../../../../imputation/cast-workflows/imputation/wdl/vntr_reference/merged_samples.sorted_clean.vcf"


#samples="samples/passing_samples_v7.1.csv"
#samples="samples/EUR_WHITE.csv"
#samples="samples/AFR_BLACK.csv"

#vcf_file=$aou_220k_imputed
#vcf_file_reheader=${vcf_file/.vcf.gz/_reheader.vcf.gz}
## Add advntr to the header, so tr_harmonizer can detect the file
#bcftools view -h $vcf_file | grep "^##" > header.txt
#echo "##source=adVNTR ver. 1.5.0" >> header.txt
#bcftools view -h $vcf_file | grep -v "^##" >> header.txt
#bcftools reheader -h header.txt $vcf_file | bcftools view -Oz > $vcf_file_reheader
#tabix -p vcf $vcf_file_reheader

# Running associatr
echo "running gwas for imputed calls"
samples="samples/AFR_BLACK.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --norm quantile \
	      --annotations annotations_acan.txt \
	      --samples $samples \
	      --norm-by-sex \
	      --is-imputed \
	      --plot
samples="samples/EUR_WHITE.csv"
./aou_gwas.py --phenotype height \
	      --num-pcs 10 \
	      --method associaTR \
	      --tr-vcf $aou_220k_imputed \
	      --norm quantile \
	      --annotations annotations_acan.txt \
	      --samples $samples \
	      --norm-by-sex \
	      --is-imputed \
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
