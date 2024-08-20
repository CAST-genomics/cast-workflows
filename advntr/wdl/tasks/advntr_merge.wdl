version 1.0

import "advntr.wdl" as advntr_one_batch

workflow run_merge_advntr {

    input {
        Array[String] per_batch_vcfs
        Array[String] per_batch_vcf_indexes
    }

    call merge_outputs {
        input:
            individual_vcfs = per_batch_vcfs,
            individual_vcf_indexes = per_batch_vcf_indexes
    }

    output {
        File merged_vcfs = merge_outputs.merged_vcfs
        File merged_vcfs_index = merge_outputs.merged_vcfs_index
    }

    meta {
        description: "This workflow merges vcf files from all batches"
    }
}

task merge_outputs {
    input {
        Array[File] individual_vcfs
        Array[File] individual_vcf_indexes
    }

    String out_prefix = "merged_samples"

    command <<<
        #mergeSTR --vcfs ~{sep=',' individual_vcfs} --out ~{out_prefix}
        bcftools merge -Oz ~{sep=',' individual_vcfs} > ~{out_prefix}.vcf.gz && tabix -p vcf ~{out_prefix}.vcf.gz
        bcftools sort -Oz ~{out_prefix}.vcf.gz > ~{out_prefix}.sorted.vcf.gz && tabix -p vcf ~{out_prefix}.sorted.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }
    output {
        File merged_vcfs = "~{out_prefix}.sorted.vcf.gz"
        File merged_vcfs_index = "~{out_prefix}.sorted.vcf.gz.tbi"
    }
}
