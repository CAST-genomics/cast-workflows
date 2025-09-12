version 1.0

workflow run_longtr {
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
        Int sleep_seconds = 0
        String? extra_longtr_args 
        Boolean separate_longtr_runs = false
    }

    call longtr {
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
          sleep_seconds=sleep_seconds,
          extra_longtr_args=extra_longtr_args,
          separate_longtr_runs=separate_longtr_runs
    }

    call sort_index_longtr{
      input :
        vcf=longtr.outfile
    }

    output {
       File outfile = sort_index_longtr.outvcf 
       File outfile_index = sort_index_longtr.outvcf_index
    }
    
    meta {
      description: "Run LongTR on multiple samples with default parameters"
    }
}

task longtr {
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
        Int sleep_seconds = 0
        Boolean separate_longtr_runs = false
        String? extra_longtr_args
    } 

    command <<<
      # If using AOU, get bamfiles from bams_file
      # Also sleep before launching LongTR
    bams_input=~{sep=',' bams}
    if [[ "~{using_batch_files}" == true ]] ; then
        bams_input=$(paste -sd, ~{bams_file})
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=~{GCS_OAUTH_TOKEN}
        sleep ~{sleep_seconds}
    fi

      # AoU crams are e.g. wgs_XXXXX.cram
    samps_flags=""
    if [[ "~{infer_samps_from_file}" == true ]] ; then
        samps=""
        for bam in $(echo ${bams_input} | sed 's/,/ /g') ; do
          samps="${samps},$(basename ${bam} .cram | sed 's/wgs_//' | cut -f1 -d'_')"
        done
        samps=${samps:1} # remove beginning excess comma
        samps_flags="--bam-samps ${samps} --bam-libs ${samps}"
    fi

    if [[ "~{separate_longtr_runs}" == false ]] ; then

    LongTR  \
        --bams ${bams_input} \
        --regions  ~{str_ref} \
        --bams  ${bams_input} \
        --fasta  ~{genome} \
        ${samps_flags} \
        --phased-bam  \
        --indel-flan-len 0 \
        --tr-vcf  ~{out_prefix}.vcf.gz

    else
        # Run longtr separately on each STR
        counter=0
        while IFS= read -r line
        do
          # Update token. expires every hour
          export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
          echo "$line" > str_${counter}.bed
          LongTR  \
            --bams ${bams_input} \
            --regions  ~{str_ref} \
            --bams  ${bams_input} \
            --fasta  ~{genome} \
            ${samps_flags} \
            --phased-bam  \
            --indel-flan-len 0 \
            --tr-vcf  ~{out_prefix}.vcf.gz
          counter=$((counter+1))
          sleep ~{sleep_seconds}
        done < ~{str_ref}
        # Concatenate the VCF files
        echo "##fileformat=VCFv4.1" > ~{out_prefix}.vcf
        for num in $(seq 0 $((counter-1)))
        do
           zcat ~{out_prefix}_${num}.vcf.gz | grep "^##command" >> ~{out_prefix}.vcf
        done
        zcat ~{out_prefix}_0.vcf.gz | grep "^#" | \
          grep -v fileformat | grep -v command >> ~{out_prefix}.vcf
        for num in $(seq 0 $((counter-1)))
        do
           zcat ~{out_prefix}_${num}.vcf.gz | grep -v "^#" >> ~{out_prefix}.vcf
        done
        bgzip ~{out_prefix}.vcf
    fi
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/hipstr-longtr:latest"
        memory: "16G"
        cpu: 2
      	maxRetries: 3
        preemptible: 3
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}

task sort_index_longtr {
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
