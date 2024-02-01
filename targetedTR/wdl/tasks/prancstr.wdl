version 1.0

workflow run_prancstr {
    input {
        File input_vcf
        String vcftype
        String readfield
        String region
        String out_prefix
    }

    call prancstr {
        input : 
          input_vcf=input_vcf,
          vcftype=vcftype,
          readfield=readfield,
          region=region,
          out_prefix=out_prefix
    }

    output {
       File outfile = prancstr.outtab
    }
    meta {
      description: "Run prancSTR on a vcf file with default parameters"
    }
}

task prancstr {
    input {
        File input_vcf
        String vcftype
        String readfield
        String region
        String out_prefix
    } 

    command <<<
      prancSTR \
          --vcf ~{input_vcf} \
          --vcftype ~{vcftype} \
          --readfield ~{readfield} \
          --out ~{out_prefix} \
          --region ~{region} \
          --only-passing 2>../../tests/prancSTR_test.log 
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-5.1.1:latest"
        memory: "8 GB"
        cpu: 1
    }

    output {
       File outtab = "${out_prefix}.tab"
    }
}
