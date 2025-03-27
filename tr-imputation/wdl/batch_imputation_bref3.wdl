version 1.0

workflow batch_imputation_bref3 {
    input {
        Array[File] snp_vcf_list = []
        Array[File] snp_idx_list = []
        File ref_panel
        File genetic_map
        Int mem
    }

    Int batch_size = length(snp_vcf_list)
    scatter(i in range(batch_size)) {
        File vcf = snp_vcf_list[i]
        File idx = snp_idx_list[i]
        String curr_file_str = basename(vcf, ".vcf.gz") 
    
        call imputation {
            input :
                snp_vcf = vcf,
                snp_idx = idx,
                ref = ref_panel,
                map = genetic_map,
                out_prefix = curr_file_str,
                mem = mem
        }
#        call extract_TR {
#            input :
#                vcf = imputation.snp_strs
#                vcf_idx = imputation.snp_strs_idx
#                out_prefix = curr_file_str
#        }
    }

    output {
        Array[File] imputed_snp_strs = imputation.snp_strs
        Array[File] imputed_snp_strs_idx = imputation.snp_strs_idx
    }
}

task imputation {
    input {
      File snp_vcf
      File snp_idx
      File ref
      File map
      String out_prefix
      Int mem
    }
    command <<<
        set -e
        ulimit -n 800000
        java -Xmx25g -jar /beagle.jar \
            gt=~{snp_vcf} \
            ref=~{ref} \
            ap=true \
            map=~{map} \
            out=~{out_prefix}
        
        tabix -p vcf ~{out_prefix}.vcf.gz 
      >>>

      runtime {
          docker: "gcr.io/ucsd-medicine-cast/beagle:latest"
          memory: mem+"G"
          maxRetries: 1
      }
      
      output {
          File snp_strs = "~{out_prefix}.vcf.gz"
          File snp_strs_idx = "~{out_prefix}.vcf.gz.tbi"
      }
}

#task extract_TR {
#    input {
#        File vcf
#        File vcf_idx
#        String out_prefix
#    }
#    
#    command <<<
#        set -e
#        ulimit -n 800000
#        echo "start extracting TRs..." 
#        bcftools view -i 'ID ~"EnsTR"' ~{vcf} -Oz -o ~{out_prefix}_imputed_TR.vcf.gz
#        echo "finish extracting, start index now" 
#        tabix -p vcf ~{out_prefix}_imputed_TR.vcf.gz
#        
#       # if [[ "~{output_snps}" == "true" ]]
#       # then
#       #     echo "start extracting SNPs..."
#       #     bcftools view -e 'ID ~"EnsTR"' ~{vcf} -Oz -o ~{out_prefix}_imputed_SNP.vcf.gz
#       #     tabix -p vcf ~{out_prefix}_imputed_SNP.vcf.gz
#       # fi 
#    >>>
#
#    output {
#        File imputed_TRs = "${out_prefix}_imputed_TR.vcf.gz"
#        File imputed_idxs = "${out_prefix}_imputed_TR.vcf.gz.tbi"
#    }
#
#}
