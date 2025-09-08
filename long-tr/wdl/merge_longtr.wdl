version 1.0

workflow merge_longtr {
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

    output {
        File outfile = mergestr.outvcf
    }

    meta {
        description: "Merge VCFs from multiple LongTR runs"
    }
}

task mergestr {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        String out_prefix
        Int total = length(vcfs)
    }

    command <<<
        touch vcf.list
        FILEARRAY=(~{sep=' ' vcfs}) # Load array into bash variable
        for (( c = 0; c < ~{total}; c++ )) # bash array are 0-indexed ;
        do
            f=${FILEARRAY[$c]}
            vcf-validator $f && echo $f >> vcf.list
            vcf-validator $f || echo "Failed: " $f
        done

        mergeSTR --vcfs-list vcf.list --vcftype longtr --out ~{out_prefix}
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-6.1.0:lates"
        memory: "4G"
    }

    output {
        File outvcf = "${out_prefix}.vcf"
    }
}
