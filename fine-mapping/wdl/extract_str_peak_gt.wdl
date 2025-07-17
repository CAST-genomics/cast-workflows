version 1.0

workflow extract_str_peak_gt {
    input {
        File pgen
        File psam
        File pvar
        File region
        String out_prefix 

    }

    call extract_peaks_str {
        input:
            pgen=pgen,
            psam=psam,
            pvar=pvar,
            region=region,
            out_prefix=out_prefix
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
        String out_prefix
    }

    String pfile_prefix = basename(pgen,".pgen")

    command <<<


        plink2 \
            --pfile ~{pfile_prefix} \
            --extract range ~{region} \
            --make-pfile \
            --out ~{out_prefix}
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"

    }

    output {
        File outfile_pgen = "${out_prefix}.pgen"
        File outfile_psam = "${out_prefix}.psam"
        File outfile_pvar = "${out_prefix}.pvar"
    }
}