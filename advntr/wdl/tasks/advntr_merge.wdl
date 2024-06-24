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
        mergeSTR --vcfs ~{sep=',' individual_vcfs} --out ~{out_prefix}
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-5.0.1:latest"
    }
    output {
        File merged_vcfs = "~{out_prefix}.vcf"
    }
}
