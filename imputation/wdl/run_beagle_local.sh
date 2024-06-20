#!/bin/bash

chr="chr15"
#gt="data/ALL.${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_bialleleic_uniq_id_name.sorted.vcf.gz"
#gt="data/ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_chr.vcf.gz"
#gt="vntr_data/acan_10mbp_aou_100_samples_no_impute_yet.vcf.gz"
#gt="vntr_data/output_chr15_acan_4mbp_aou_275_lrwgs_samples.vcf.gz"
#gt="vntr_data/acan_10mbp_aou_100_samples_no_impute_yet.vcf.gz"
#gt="vntr_data/acan_10mbp_aou_185_lrwgs_samples_subset.vcf.gz"
gt="vntr_data/output_chr15_acan_10mbp_aou_10k_samples.vcf.gz"
#ref="data/${chr}_final_SNP_merged_additional_TRs.vcf.gz"
#ref="data/${chr}_final_SNP_merged_additional_TRs.bref3"
#ref="vntr_data/phased_ACAN_vntr_snp_650_samples.sorted.vcf.gz"
#ref="../../../../imputation/cast-workflows/imputation/wdl/vntr_reference/phased_ACAN_vntr_snp_740_samples.sorted.vcf.gz"
ref="../../../../imputation/cast-workflows/imputation/wdl/vntr_reference/ref_phased_reheader_chr15_acan_vntr.vcf.gz"

out="data/output_${chr}_10mb_aou_10k_srwgs_w8_o2"
#out="data/output_${chr}_4mb_aou_275_lrwgs"

beagle_5_4_old="beagle.19Apr22.7c0.jar"
beagle=$beagle_5_4_old
java -Xmx4g -jar $beagle \
        gt=$gt \
        ref=$ref \
        out=$out \
        window=8 \
        overlap=2 \
	chrom=$chr && \
tabix -p vcf $out.vcf.gz
