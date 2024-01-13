version 1.0

workflow run_hipstr {
    input {
        Array[File] bams = []
        Array[File] bam_indices = []
        File genome
        File genome_index
        File str_ref
        String out_prefix
        Boolean infer_samps_from_file = false
        Boolean using_batch_files = false
        File? bams_file
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        Int? sleep_seconds = 0
    }

    call hipstr {
        input : 
          bams=bams, 
          bam_indices=bam_indices,
          genome=genome, 
          genome_index=genome_index,
          str_ref=str_ref,
          out_prefix=out_prefix,
          infer_samps_from_file=infer_samps_from_file,
          using_batch_files=using_batch_files,
          bams_file=bams_file,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
          sleep_seconds=sleep_seconds
    }

    call sort_index_hipstr {
      input :
        vcf=hipstr.outfile
    }

    output {
       File outfile = sort_index_hipstr.outvcf 
       File outfile_index = sort_index_hipstr.outvcf_index
    }
    
    meta {
      description: "Run HipSTR on multiple samples with default parameters"
    }
}

task hipstr {
    input {
        Array[File] bams
        Array[File] bam_indices
        File genome
        File genome_index
        File str_ref
        String out_prefix
        Boolean infer_samps_from_file = false
        Boolean using_batch_files = false
        File? bams_file
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        Int? sleep_seconds = 0
    } 

    command <<<
      samps_flags=""
      if [[ "~{infer_samps_from_file}" == true ]] ; then
        samps=""
        for bam in ~{sep=" " bams} ; do
          samps="${samps},$(basename "${bam}" | cut -f 1 -d _)"
        done
        samps=${samps:1} # remove beginning excess comma
        samps_flags="--bam-samps ${samps} --bam-libs ${samps}"
      fi

      # If using AOU, get bamfiles from bams_file
      # Also sleep before launching HipSTR
      bams_input=~{sep=',' bams}
      if [[ "~{using_batch_files}" == true ]] ; then
        bams_input=$(paste -sd, ~{bams_file})
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=~{GCS_OAUTH_TOKEN}
        sleep ${sleep_seconds}
      fi
      HipSTR \
          --bams ${bams_input} \
          --fasta ~{genome} \
          --regions ~{str_ref} \
          --str-vcf ~{out_prefix}.vcf.gz \
          --min-reads 10 \
          ${samps_flags}
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/hipstr-gymreklab-gcs"
        memory: "16 GB"
        memoryMin: "16 GB"
        cpu: 2
      	maxRetries: 3
        preemptible: 3
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}

task sort_index_hipstr {
  input {
    File vcf
  }

  String basename = basename(vcf, ".vcf.gz")

  command <<<
    zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
  >>>

  runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

  output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}