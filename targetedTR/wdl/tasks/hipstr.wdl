version 1.0

workflow run_hipstr {
    input {
        File bam
        File bam_index
        File genome
        File genome_index
        File str_ref
        String out_prefix
    }

    call hipstr {
        input : 
          bam=bam, 
          bam_index=bam_index,
          genome=genome, 
          genome_index=genome_index,
          str_ref=str_ref,
          out_prefix=out_prefix
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
      description: "Run HipSTR on a single sample with default parameters"
    }
}

task hipstr {
    input {
        File bam
        File bam_index
        File genome
        File genome_index
        File str_ref
        String out_prefix
    } 

    command <<<
      HipSTR \
          --bams ~{bam} \
          --fasta ~{genome} \
          --regions ~{str_ref} \
          --str-vcf ~{out_prefix}.vcf.gz \
          --def-stutter-model \
          --min-reads 10
    >>>
    
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/hipstr-gymreklab:latest"
	memory: "16 GB"
        cpu: "1"
        disks: "local-disk 30 HDD"
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
