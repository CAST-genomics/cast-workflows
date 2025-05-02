version 1.0

workflow bgenTopgen {
    input {
        File bgen_file
        File sample_file
        Int split_mem 
        String chrom
        String bgen_qc
    }

    call convert {
        input :
            input_file = bgen_file,
            input_sample = sample_file,
            outprefix = chrom,
            mem = split_mem,
            qc_option = bgen_qc
    }

    output {
        Array[File] pgen_files = convert.pgen_files
    }
}

task convert {

    input {
        File input_file
        File input_sample
        String outprefix
        Int mem
        String qc_option
    }

    command <<<
        set -e
        ulimit -n 800000 
        echo "start conversion"
        
        mkdir -p bgenTopgen
        plink2 --bgen ~{input_file} "ref-first"\
            --sample ~{input_sample} ~{qc_option} \
            --make-pgen \
            --out "bgenTopgen/~{outprefix}_filtered"
        echo "all done!" 
        ls bgenTopgen
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem + "GB" 
        maxRetries: 1 
    }

    output {
        Array[File] pgen_files = glob("bgenTopgen/~{outprefix}_filtered.p*")
    }
}
