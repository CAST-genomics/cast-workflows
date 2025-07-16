version 1.0

workflow find_peaks {
    input {
        Array[File] pgens = [] 
        Array[File] psams = [] 
        Array[File] pvars = []  
        File region
        String out_prefix 

    }

    call find_peaks_str {
        input:
            pgens=pgens,
            psams=psams,
            pvars=pvars,
            region=region,
            out_prefix=out_prefix
    }
    output {
        Array[Array[File]] outfile_pgen = find_peaks_str.pgen
        Array[Array[File]] outfile_psam = find_peaks_str.psam
        Array[Array[File]] outfile_pvar = find_peaks_str.pvar
    }

    meta {
        description: "This workflow is to extract fine-mapping peak region from STR imputation STR pfile"

    }
}

task find_peaks_str {
    input {
        Array[File] pgens = [] 
        Array[File] psams = [] 
        Array[File] pvars = []  
        File region
        String out_prefix
    }
    command <<<
        PFILEARRAY=(~{sep=" " pgens})
        outfiles=""
        # bash array are 0-indexed 
        for (( c = 0; c < ~{total}; c++ ))
        do
            pfile=${PFILEARRAY[$c]}
            pfile_outprefix="${pfile%.pgen}"
            chrom_outprefix=$(basename $pfile .pgen)
            plink2 \
                --pfile ~{pfile} \
                --extract range ~{region} \
                --make-pfile \
                --out "~{out_prefix}_${chrom_outprefix}"
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        memory: "6G"
        disks: "local-disk 400 SSD"
    }

    output {
        Array[Array[File]] pgen = "${out_prefix}_${chrom_outprefix}.pgen"
        Array[Array[File]] psam = "${out_prefix}_${chrom_outprefix}.psam"
        Array[Array[File]] pvar = "${out_prefix}_${chrom_outprefix}.pvar"
    }
}