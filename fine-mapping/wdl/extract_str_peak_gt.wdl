version 1.0

workflow extract_str_peak_gt {
    input {
        File pgen
        File psam
        File pvar
        File region
        String pheno
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""

    }

    call extract_peaks_str {
        input:
            pgen=pgen,
            psam=psam,
            pvar=pvar,
            region=region,
            pheno=pheno,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN
    }
    output {
        File outfile_pgen = extract_peaks_str.outfile_pgen
        File outfile_psam = extract_peaks_str.outfile_psam
        File outfile_pvar = extract_peaks_str.outfile_pvar
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
        File region
        String pheno
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
    }


    command <<<
        pgen_outfile = ""
        pvar_outfile = ""
        psam_outfile = ""

        pfile_outprefix="${pgen%.pgen}"
        while IFS=$' ' read -r chr from to
        do
            plink2 --pfile "${pfile_prefix}" \
                   --chr "$chr" \
                   --from-bp "$from" \
                   --to-bp "$to" \
                   --make-pgen \
                   --out "$chr_$from_$to_~{pheno}"

            pgen_outfile+="$chr_$from_$to_~{pheno}.pgen "
            pvar_outfile+="$chr_$from_$to_~{pheno}.pvar "
            psam_outfile+="$chr_$from_$to_~{pheno}.psam "
                
            
        done < ${region}

        ls ${pgen_outfile}|sed -E 's/\.pgen.*//' | sort -u > "~{pheno}_$chr_str.txt"

        plink2 --pmerge-list "~{pheno}_$chr_str.txt" --sort-vars --make-pgen --out "~{pheno}_$chr_str_merged"


    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"

    }

    output {
        File outfile_pgen = "${pheno}_$chr_str_merged.pgen"
        File outfile_psam = "${pheno}_$chr_str_merged.psam"
        File outfile_pvar = "${pheno}_$chr_str_merged.pvar"
    }
}