version 1.0

workflow rerun_annotator {
    input {
        Array[File] vcf = []
        Array[File] vcf_index = []
        Array[File] ref_vcf = []
        Array[File] ref_index = []
        String out_prefix 

    }
    ### Separate workflow for each chrom vcf ###

    Int num_vcf = length(vcf)

    scatter (i in range(num_vcf)) {
        File chrom_vcf = vcf[i]
        File chrom_vcf_index = vcf_index[i]
        File chrom_ref = ref_vcf[i]
        File chrom_ref_index = ref_index[i]
        String chrom = basename(chrom_vcf,"_TR_merged.vcf.gz")

        call annotaTR {
            input:
                vcf=chrom_vcf,
                vcf_index=chrom_vcf_index,
                ref_vcf=chrom_ref,
                ref_index=chrom_ref_index,
                out_prefix=chrom+out_prefix
        }
    }
    output {
        Array[File] outfile_pgen = annotaTR.pgen
        Array[File] outfile_psam = annotaTR.psam
        Array[File] outfile_pvar = annotaTR.pvar
        Array[File] outfile_vcf = annotaTR.outvcf
        Array[File] outfile_vcfind = annotaTR.outvcfind
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
            --outtype pgen vcf \
            --vcf-outtype z \
            --dosages beagleap_norm \
            --ignore-duplicates \
            --match-refpanel-on locid \
            --warn-on-AP-error \
            --update-ref-alt

        tabix -p vcf ~{out_prefix}.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-annotatr:yli091230"
        disks: "local-disk 200 SSD"
        memory: "40G"
    }

    output {
        File pgen = "${out_prefix}.pgen"
        File psam = "${out_prefix}.psam"
        File pvar = "${out_prefix}.pvar"
        File outvcf = "${out_prefix}.vcf.gz"
        File outvcfind = "${out_prefix}.vcf.gz.tbi"
    }
}