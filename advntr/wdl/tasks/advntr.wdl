version 1.0

import "advntr_single.wdl" as advntr_single

workflow run_advntr {

    input {
        Array[String] bam_files
        String region_file
        String google_project
        String vntr_id
        Int mem
    }

    scatter (i in range(length(bam_files))) {
        String bam_file = bam_files[i]
        call advntr_single.advntr_single_sample as advntr_single_sample {
            input:
                bam_file=bam_file,
                region_file=region_file,
                google_project=google_project,
                vntr_id=vntr_id,
                sleep_seconds=i,
                mem=mem
        }
    }

    call merge_sort {
        input:
            individual_vcfs = advntr_single_sample.out_vcf,
            individual_vcf_indexes = advntr_single_sample.out_vcf_index,
            mem = mem*2
    }

    output {
        File merged_vcf = merge_sort.merged_vcf
        File merged_vcf_index = merge_sort.merged_vcf_index
    }

    meta {
        description: "This workflow calls adVNTR to genotype VNTRs for a single batch in parallel"
    }
}

task merge_sort {
    input {
        Array[File] individual_vcfs
        Array[File] individual_vcf_indexes
        Int mem
    }

    String out_prefix = "merged_samples"

    command <<<
        touch ~{sep=' ' individual_vcf_indexes}
        echo "Calling merge"
        date
        bcftools merge --force-single --merge id -Oz ~{sep=' ' individual_vcfs} > ~{out_prefix}.vcf.gz
        date
        echo "Merged"
        echo "Indexing"
        date
        tabix -p vcf ~{out_prefix}.vcf.gz
        date
        echo "Indexed"
        echo "Calling Sort"
        date
        bcftools sort -Oz ~{out_prefix}.vcf.gz > ~{out_prefix}.sorted.vcf.gz && tabix -p vcf ~{out_prefix}.sorted.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        disks: "local-disk ${mem} SSD"
    }

    output {
        File merged_vcf = "~{out_prefix}.sorted.vcf.gz"
        File merged_vcf_index = "~{out_prefix}.sorted.vcf.gz.tbi"
    }
}
