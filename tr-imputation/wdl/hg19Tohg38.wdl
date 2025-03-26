version 1.0

workflow hg19Tohg38 {
    input {
        File hg19_vcf
        Int liftover_mem
        }
    call liftover {
        input :
            input_vcf = hg19_vcf,
            mem = liftover_mem 
    } 
    output {
        File hg38_vcf = liftover.hg38_vcfs
        File hg38_idx = liftover.hg38_idxs
        File unlifted_vcf = liftover.unlifted_vcfs
        File logs = liftover.liftover_logs
    }
}

task liftover {
      input {
         File input_vcf
         Int mem
      }
      String vcf_prefix = basename(input_vcf, ".vcf.gz") 
      command <<<
          set -e
          echo "start liftover ~{vcf_prefix}"

          liftover.py "~{input_vcf}" .
          ls -lht
          mkdir ./sort_temp
          max_mem=$((~{mem}-4))
          bcftools sort -Oz -T ./sort_temp -m "${max_mem}G" -o ./~{vcf_prefix}_hg38_sorted.vcf.gz ./~{vcf_prefix}_hg38.vcf.gz
          tabix -p vcf ./~{vcf_prefix}_hg38_sorted.vcf.gz
      >>>

      runtime {
          docker: "gcr.io/ucsd-medicine-cast/pyliftover:latest"
          memory: mem + "GB"
          disk: "local-disk 200 SSD"
      }

      output {
          File hg38_vcfs="./${vcf_prefix}_hg38_sorted.vcf.gz"
          File hg38_idxs="./${vcf_prefix}_hg38_sorted.vcf.gz.tbi"
          File unlifted_vcfs="./${vcf_prefix}_unlifted.vcf.gz"
          File liftover_logs="./${vcf_prefix}_liftOver.log"
      }

}
