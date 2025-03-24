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
        Array[File] hg38_vcf = liftover.hg38_vcfs
        Array[File] hg38_idx = liftover.hg38_idxs
        Array[File] unlifted_vcf = liftover.unlifted_vcfs
        Array[File] logs = liftover.liftover_logs
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
          tabix -p vcf ./~{vcf_prefix}_hg38.vcf.gz
      >>>

      runtime {
          docker: "gcr.io/ucsd-medicine-cast/pyliftover:latest"
          memory: mem + "GB"
          disk: "local-disk 200 SSD"
      }

      output {
          Array[File] hg38_vcfs=glob("./*hg38.vcf.gz")
          Array[File] hg38_idxs=glob("./*hg38.vcf.gz.tbi")
          Array[File] unlifted_vcfs=glob("./*unlifted.vcf.gz")
          Array[File] liftover_logs=glob("./*liftOver.log")
      }

}
