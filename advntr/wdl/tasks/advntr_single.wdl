
version 1.0

workflow advntr_single_sample {
    input {
        String bam_file
        String region
        String google_project
        String gcloud_token
        String vntr_id
        Int sleep_seconds
    }
    call download_input {
        input :
        bam_file = bam_file,
        region = region,
        google_project = google_project,
        gcloud_token = gcloud_token,
    }
    call genotype {
        input :
        region = region,
        vntr_id = vntr_id,
        target_bam_file = download_input.target_bam_file,
        target_bam_index_file = download_input.target_bam_index_file,
        sleep_seconds = sleep_seconds
    }
    
    call sort_index {
        input :
        vcf = genotype.genotype_output
    }
    output {
        File out_vcf = sort_index.out_vcf
        File out_vcf_index = sort_index.out_vcf_index
    }
    meta {
        description: "This workflow calls adVNTR to genotype VNTRs for a single sample"
    }
}



task sort_index {
  input {
    File vcf
  }

  String basename = basename(vcf, ".vcf")

  command <<<
    cat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
  >>>

  runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
         maxRetries: 3
    }

  output {
    File out_vcf = "${basename}.sorted.vcf.gz"
    File out_vcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}

task download_input {
    input {
        String bam_file
        String region
        String google_project
        String gcloud_token
    }


    # Names of all the intermediate files generated to get the target bam file.
    String unsorted_target_bam = "target_region_~{sample_id}_unsorted.bam"
    String sorted_target_bam = "target_~{sample_id}.bam"
    String sorted_target_bam_index = "target_~{sample_id}.bam.bai"
    String sample_id = sub(basename(bam_file), ".bam", "")

    #/google-cloud-sdk/bin/gcloud init --skip-diagnostics
    command <<<
        ls -lh .
        export HTSLIB_CONFIGURE_OPTIONS="--enable-gcs"
        echo "pwd $(pwd)"
        /google-cloud-sdk/bin/gcloud --version
        /google-cloud-sdk/bin/gcloud config list --format='text(core.project)'
        export gcloud_token=$(/google-cloud-sdk/bin/gcloud auth application-default print-access-token --project ~{google_project})
        export GCS_OAUTH_TOKEN=${gcloud_token}
        export GCS_REQUESTER_PAYS_PROJECT="~{google_project}"
        samtools view -hb -o ~{unsorted_target_bam} --use-index ~{bam_file} ~{region}
        samtools sort -o ~{sorted_target_bam} ~{unsorted_target_bam}
        samtools index ~{sorted_target_bam}
        ls -lh .
    >>>

    runtime {
        docker:"sarajava/samtools:1.13_gcli"
        maxRetries: 3
        memory: "4G"
    }
    output {
        File target_bam_file = "~{sorted_target_bam}"
        File target_bam_index_file = "~{sorted_target_bam_index}"
    }
}

task genotype {
    input {
        String region
        String vntr_id
        File target_bam_file
        File target_bam_index_file
        Int sleep_seconds
    }

    # Provide the names of all the output files being generated by AdVNTR including intermediate files.
    String logging = "./log_~{sample_id}.bam.log"
    String filtering_out = "./filtering_out_~{sample_id}.unmapped.fasta.txt"
    String keywords = "./keywords_~{sample_id}.unmapped.fasta.txt"
    String unmapped = "./~{sample_id}.unmapped.fasta"
    String genotype_output = "./~{sample_id}.vcf"
    String sample_id = sub(basename(target_bam_file), ".bam", "")


    # VNTR_db is placed in the docker file. So the path is within the docker image.
    String vntr_db = "/adVNTR/vntr_db/hg38_VNTRs_by_TRF.db"
    #String reference_name = "/adVNTR/hg38_reference/hg38full.fa"

    # Set VNTR ids for genotyping based on input.
    # Two options right now: VNTR in the ACAN gene or the list of 52 phenotype associated VNTRs.
    #String vids = "$(cat /adVNTR/vntr_db/phenotype_associated_vntrs_comma.txt)"

    command <<<
        sleep ~{sleep_seconds}
        echo "num reads $(samtools view -c ~{target_bam_file})"
        echo "num reads in region $(samtools view -c ~{target_bam_file} ~{region}) "
        /usr/bin/time -v advntr genotype  \
        --alignment_file ~{target_bam_file} \
        --models ~{vntr_db}  \
        --working_directory . \
        -vid ~{vntr_id} \
        --outfmt vcf \
        --log_pacbio_reads \
        --pacbio > ~{genotype_output}
    >>>

    runtime {
        docker:"sarajava/advntr:1.5.0_v12"
        memory: "4G"
        maxRetries: 3
    }

    output {
        File? log_file = "~{logging}"
        File? filtering_out = "~{filtering_out}"
        File? keywords = "~{keywords}"
        File genotype_output = "~{genotype_output}"
    }
}
