version 1.0

workflow rerun_annotator {
    input {
        File vcf 
        File vcf_index 
        File ref_vcf 
        File ref_index 
        String out_prefix 

    }

    call annotaTR {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            ref_vcf=ref_vcf,
            ref_index=ref_index,
            out_prefix=out_prefix+"_reannotated"
    }

    output {
        File outfile_pgen = annotaTR.pgen
        File outfile_psam = annotaTR.psam
        File outfile_pvar = annotaTR.pvar
    }

    meta {
        description: "This workflow reruns annotaTR on imputation TR outputs to output a new INFO filed named DSCOUNT contains the non-zero dosage and their counts number in the pvar file."
    }                 
}

task annotaTR {
    input {
        File vcf 
        File vcf_index
        File ref_vcf
        File ref_index
        String out_prefix
    }
    
    command <<<
        set -e
        annotaTR --vcf ~{vcf} \
            --ref-panel ~{ref_vcf} \
            --out ~{out_prefix} \
            --vcftype hipstr \
            --outtype pgen \
            --vcf-outtype z \
            --dosages beagleap_norm \
            --ignore-duplicates \
            --match-refpanel-on locid \
            --warn-on-AP-error \
            --update-ref-alt

    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-annotatr:yli091230"
        disks: "local-disk 120 SSD"
        memory: "30G"
    }

    output {
        File pgen = "${out_prefix}.pgen"
        File psam = "${out_prefix}.psam"
        File pvar = "${out_prefix}.pvar"
    }
}