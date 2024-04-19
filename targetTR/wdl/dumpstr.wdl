version 1.0

workflow run_dumpstr {
    input {
        File vcf
        String out_prefix
    }

    call dumpstr {
        input : 
          vcf=vcf,
          out_prefix=out_prefix
    }

    output {
       File outfile = dumpstr.outvcf 
    }
    
    meta {
      description: "Run dumpSTR on a (HipSTR) VCF"
    }
}

task dumpstr {
  input {
    File vcf
    String out_prefix
  }

  command <<<
      dumpSTR --vcf ~{vcf} --out ~{out_prefix} \
        --hipstr-min-call-Q 0.9 \
        --hipstr-min-call-DP 10 \
        --hipstr-max-call-DP 10000 \
        --hipstr-min-supp-reads 2 \
        --hipstr-max-call-stutter 0.15 \
        --hipstr-max-call-flank-indel 0.15
  >>>
    
  runtime {
      docker: "gcr.io/ucsd-medicine-cast/trtools-6.0.1:latest"
  }

  output {
      File outvcf = "${out_prefix}.vcf"
  }
}
