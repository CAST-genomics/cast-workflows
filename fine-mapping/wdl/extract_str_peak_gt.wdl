version 1.0

workflow extract_str_peak_gt {
    input {
        File pgen
        File psam
        File pvar
        String pheno
        File tr_list
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""

    }

    call extract_peaks_str {
        input:
            pgen=pgen,
            psam=psam,
            pvar=pvar,
            pheno=pheno,
            tr_list=tr_list,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN
    }
    output {
        Array[Array[File]] outfile_pgen = extract_peaks_str.outfile_pgen
        Array[Array[File]] outfile_psam = extract_peaks_str.outfile_psam
        Array[Array[File]] outfile_pvar = extract_peaks_str.outfile_pvar
    }

    meta {
        description: "This workflow is to extract fine-mapping peak region from STR imputation pfile"

    }
}

task extract_peaks_str {
    input {
        File pgen
        File psam
        File pvar
        String pheno
        File tr_list
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
    }


    command <<<
        set -e
        pfile_outprefix="${pgen%.pgen}"

        while IFS=$' ' read -r -a locus
        do
            chr=${locus[0]}
            from=${locus[1]}
            to=${locus[2]}
            
            plink2 --pfile "${pfile_outprefix}"  \
                --chr "$chr" \
                --from-bp "$from" \
                --to-bp "$to" \
                --make-pgen \
                --out "${pheno}_${chr}_${from}_${to}_str"     
            
        done < ~{tr_list}

    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"

    }

    output {
        Array[Array[File]] outfile_pgen = "${pheno}_$chr_$from_$to_str.pgen"
        Array[Array[File]] outfile_psam = "${pheno}_$chr_$from_$to_str.psam"
        Array[Array[File]] outfile_pvar = "${pheno}_$chr_$from_$to_str.pvar"
    }
}