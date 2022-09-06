version 1.0

workflow run_HipSTR_single {
    input {
        File bam = "/Users/yang/Desktop/CAST_repository/cast-workflows/hipstr/wdl/tasks/tests/HG000190_CBL_sorted.cram"
        File genome = "/Users/yang/Desktop/CAST_repository/cast-workflows/hipstr/wdl/tasks/tests/chr11.fa"
        File str_ref= "/Users/yang/Desktop/CAST_repository/cast-workflows/hipstr/wdl/tasks/tests/hipref_test.bed"
        File bam_index = "/Users/yang/Desktop/CAST_repository/cast-workflows/hipstr/wdl/tasks/tests/HG000190_CBL_sorted.crai"
        File stutter_model
    }

    call hipstr {
        input : 
          bam=bam, 
          genome=genome, 
          str_ref=str_ref,
          bam_index=bam_index,
          stutter_model=stutter_model
    }

    output {
       File outfile = hipstr.outfile 
    }
    meta {
      description: "This workflow use HipSTR call STRs from single sample based on pre-generated stutter model file"
    }
}

task hipstr {
    input {
        File bam
        File genome
        File str_ref
        File bam_index
        File stutter_model
    } 

    String vcffile=basename(bam)
    command <<<
      echo ~{vcffile} 
      HipSTR \
          --bams ~{bam} \
          --fasta ~{genome} \
          --regions ~{str_ref} \
          --str-vcf ~{vcffile}.vcf.gz \
          --stutter-in ~{stutter_model}
    >>>
    
    runtime {
        docker:"yli091230/hipstr:small"
    }

    output {
       File outfile = "${vcffile}.vcf.gz"
    }

}



