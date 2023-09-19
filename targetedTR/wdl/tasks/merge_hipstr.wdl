version 1.0

workflow merge_hipstr {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        String out_prefix
    }

    call mergestr {
        input : 
          vcfs=vcfs,
          vcf_indexes=vcf_indexes,
          out_prefix=out_prefix+".merged"
    }
# call mergefix only when there are multiple records for one repeat and change outfile in the output, cant have two output have the same name
    #call mergefix {
    #    input :
    #      vcf=mergestr.outvcf,
    #      outfile=out_prefix+".merged.vcf"
    #}

    output {
       #File outfile = mergefix.outvcf 
       File outfile = mergestr.outvcf
    }

    meta {
      description: "Merge VCFs from multiple HipSTR runs"
    }
}

task mergestr {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String out_prefix
  }

  command <<<
      mergeSTR --vcfs ~{sep=',' vcfs} --out ~{out_prefix}
  >>>
    
  runtime {
      docker: "gcr.io/ucsd-medicine-cast/trtools-5.0.1:latest"
  }

  output {
      File outvcf = "${out_prefix}.vcf"
  }
}

task mergefix {
  input {
    File vcf
    String outfile
  }

  command <<<
    Hipstr_correction.py ~{vcf} ~{outfile}
  >>>

  runtime {
    docker: "gcr.io/ucsd-medicine-cast/ensembletr:latest"
  }

  output {
    File outvcf = "${outfile}"
  }
}
