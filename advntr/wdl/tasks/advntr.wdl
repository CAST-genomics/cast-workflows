version 1.0

import "advntr_single.wdl" as advntr_single

workflow run_advntr {

    input {
        Array[String] bam_files
        String region
        String google_project
        String gcloud_token
        String vntr_id
    }

    scatter (i in range(length(bam_files))) {
        String bam_file = bam_files[i]
        call advntr_single.advntr_single_sample as advntr_single_sample {
            input:
                bam_file=bam_file,
                region=region,
                google_project=google_project,
                gcloud_token=gcloud_token,
                vntr_id=vntr_id,
                sleep_seconds=i,
        }
    }

    call merge_outputs {
        input:
            individual_vcfs = advntr_single_sample.out_vcf,
            individual_vcf_indexes = advntr_single_sample.out_vcf_index
    }

    call advntr_single.sort_index as sort_index {
        input:
            vcf=merge_outputs.merged_vcfs
    }

    output {
        File merged_vcf = sort_index.out_vcf
        File merged_vcf_index = sort_index.out_vcf_index
    }

    meta {
        description: "This workflow calls adVNTR to genotype VNTRs for a single batch in parallel"
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
