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
        Int sleep_seconds = 0
        String? extra_hipstr_args = "--min-reads 10"
        Boolean separate_hipstr_runs = false
        Boolean longtr = false
	      Int hipstr_mem = 16

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
          sleep_seconds=sleep_seconds,
          hipstr_mem=hipstr_mem,
          extra_hipstr_args=extra_hipstr_args,
          longtr=longtr,
          separate_hipstr_runs=separate_hipstr_runs
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
        Int sleep_seconds = 0
        Boolean separate_hipstr_runs = false
        Boolean longtr = false
	      Int hipstr_mem = 16
        String? extra_hipstr_args = "--min-reads 10"
       
    } 

    command <<<
      # If using AOU, get bamfiles from bams_file
      # Also sleep before launching HipSTR
      bams_input=~{sep=',' bams}
      if [[ "~{using_batch_files}" == true ]] ; then
        bams_input=$(paste -sd, ~{bams_file})
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=~{GCS_OAUTH_TOKEN}
        sleep ~{sleep_seconds}
      fi

      # AoU crams are e.g. wgs_XXXXX.cram
      # UKB crams are ???
      samps_flags=""
      if [[ "~{infer_samps_from_file}" == true ]] ; then
        samps=""
        for bam in $(echo ${bams_input} | sed 's/,/ /g') ; do
          samps="${samps},$(basename ${bam} .cram | sed 's/wgs_//' | cut -f1 -d'_')"
        done
        samps=${samps:1} # remove beginning excess comma
        samps_flags="--bam-samps ${samps} --bam-libs ${samps}"
      fi

      if [[ "~{longtr}" == false ]] ; then
        if [[ "~{separate_hipstr_runs}" == false ]] ; then
        HipSTR \
            --bams ${bams_input} \
            --fasta ~{genome} \
            --regions ~{str_ref} \
            --str-vcf ~{out_prefix}.vcf.gz \
            ${samps_flags} \
            ~{extra_hipstr_args} 
        else
          # Run HipSTR separately on each STR
          counter=0
          while IFS= read -r line
          do
            # Update token. expires every hour
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            echo "$line" > str_${counter}.bed
            HipSTR \
              --bams ${bams_input} \
              --fasta ~{genome} \
              --regions str_${counter}.bed \
              --str-vcf ~{out_prefix}_${counter}.vcf.gz \
              ${samps_flags} \
              ~{extra_hipstr_args} 
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
      
      else
        LongTR --help

          #LongTR \
          #      --bams  ${bams_input} \         
          #      --fasta  ~{genome} \
          #      --regions  ~{str_ref} \
          #      --tr-vcf  ~{out_prefix}.vcf.gz \
          #      --phased-bam \
          #      ${samps_flags} 

  
      fi
     
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/hipstr-longtr:latest"
        memory: hipstr_mem + "GB"
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
