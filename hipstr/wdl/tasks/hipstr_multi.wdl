version 1.0

workflow run_HipSTR_multi {
    input {
        Array[File] bams 
        File genome 
        File str_ref
        Array[File] bam_index 
    }

    call hipstr_multi {
        input : 
          bams=bams, 
          genome=genome, 
          str_ref=str_ref,
          bam_index=bam_index
    }

    output {
       File outfile = hipstr_multi.outfile
       File stutter_file = hipstr_multi.stutter_file
    }
    meta {
      description: "This workflow use HipSTR call STRs from multiple sample and output the stutter model file for running single sample"
    }
}

task hipstr_multi {
    input {
        Array[File] bams
        File genome
        File str_ref
        Array[File] bam_index
    } 

    String vcffile="joint_call"
    command <<<
      echo ~{vcffile} 
      HipSTR \
          --bams ~{sep=',' bams} \
          --fasta ~{genome} \
          --regions ~{str_ref} \
          --str-vcf ~{vcffile}.vcf.gz \
          --stutter-out stutter_models.txt
    >>>
    
    runtime {
        docker:"yli091230/hipstr:small"
    }

    output {
       File outfile = "${vcffile}.vcf.gz"
       File stutter_file = "stutter_models.txt"
    }

}



