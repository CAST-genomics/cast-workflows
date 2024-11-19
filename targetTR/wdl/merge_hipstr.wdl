version 1.0

workflow merge_hipstr {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        String out_prefix
	      Int? merge_mem = 4
        Boolean longtr = false
    }

    call mergestr {
        input : 
          vcfs=vcfs,
          vcf_indexes=vcf_indexes,
          out_prefix=out_prefix+".merged",
	        merge_mem=merge_mem,
          longtr=longtr
    }

    output {
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
    Int total = length(vcfs)
    Int merge_mem = 4
    Boolean longtr = false
  }

  command <<<
      touch vcf.list
      FILEARRAY=(~{sep=' ' vcfs}) # Load array into bash variable
      for (( c = 0; c < ~{total}; c++ )) # bash array are 0-indexed ;
      do
           f=${FILEARRAY[$c]}
           vcf-validator $f && echo $f >> vcf.list
           vcf-validator $f || echo "Failed: " $f
      #done
      echo $vcf.list
      if [[ "~{longtr}" == false ]] ; then
        mergeSTR --vcfs-list vcf.list --out ~{out_prefix}
      else
        mergeSTR --vcftype longtr --vcfs-list vcf.list --out ~{out_prefix}
  >>>
    
  runtime {
      #docker: "gcr.io/ucsd-medicine-cast/trtools-mergestr-files:latest"
      docker: "gcr.io/ucsd-medicine-cast/trtools-6.0.2:latest"
      memory: merge_mem +"GB"
  }

  output {
      File outvcf = "${out_prefix}.vcf"
  }
}
