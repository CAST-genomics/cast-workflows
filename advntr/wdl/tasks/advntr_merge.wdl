version 1.0

workflow run_merge_advntr {

    input {
        Array[String] per_batch_vcfs
        Array[String] per_batch_vcf_indexes
        Int mem
    }

    call merge_outputs {
        input:
            individual_vcfs = per_batch_vcfs,
            individual_vcf_indexes = per_batch_vcf_indexes,
            mem=mem
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
        Int mem
    }

    String out_prefix = "merged_samples_batches"

    command <<<
        touch ~{sep=' ' individual_vcf_indexes}
        echo "Merging vcfs"
        bcftools merge --force-samples --merge id -Oz ~{sep=' ' individual_vcfs} > ~{out_prefix}.vcf.gz && tabix -p vcf ~{out_prefix}.vcf.gz
        echo "Sorting and indexing vcfs"
        bcftools sort -Oz ~{out_prefix}.vcf.gz > ~{out_prefix}.sorted.vcf.gz && tabix -p vcf ~{out_prefix}.sorted.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        disks: "local-disk ${mem} SSD"
    }
    output {
        File merged_vcfs = "~{out_prefix}.sorted.vcf.gz"
        File merged_vcfs_index = "~{out_prefix}.sorted.vcf.gz.tbi"
    }
}
