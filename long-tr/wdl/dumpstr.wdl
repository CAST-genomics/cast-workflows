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
      description: "Run dumpSTR on a (LongTR) VCF"
    }
}

task dumpstr {
    input {
        File vcf
        String out_prefix
    }

    command <<<
            dumpSTR --vcf ~{vcf} --out ~{out_prefix} \
                --longtr-min-call-Q 0.9 \
                --longtr-min-call-DP 10 \
                --longtr-max-call-DP 10000 \
                --longtr-min-supp-reads 2 \
                --longtr-max-call-flank-indel 0.15
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-6.1.0:latest"
    }

    output {
        File outvcf = "${out_prefix}.vcf"
    }
}
