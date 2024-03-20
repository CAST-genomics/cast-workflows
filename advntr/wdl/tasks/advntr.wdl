version 1.0

import "run_advntr_single.wdl" as run_advntr_single

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
        Int sleep_seconds = i
        call run_advntr_single.advntr_single_sample as advntr_single_sample {
            input:
                bam_file=bam_file,
                region=region,
                google_project=google_project,
                gcloud_token=gcloud_token,
                vntr_id=vntr_id,
                sleep_seconds=sleep_seconds,
        }
    }

    call merge_outputs {
        input:
            individual_vcfs = advntr_single_sample.out_vcf,
            individual_vcf_indexes = advntr_single_sample.out_vcf_index
    }

    output {
        File merged_vcfs = merge_outputs.merged_vcfs
    }

    meta {
        description: "This workflow calls adVNTR to genotype VNTRs"
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
        #cpu: "1"
        #memory: "4G"
    }
    output {
        File merged_vcfs = "~{out_prefix}.vcf"
    }
}

