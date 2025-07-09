version 1.0

workflow convert_vcf_to_pgen {
    input {
        File vcf 
        File vcf_index 
        String out_prefix 

    }

    call convert_vcf_to_pgen {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            out_prefix=out_prefix
    }
    output {
        File outfile_pgen = convert_vcf_to_pgen.pgen
        File outfile_psam = convert_vcf_to_pgen.psam
        File outfile_pvar = convert_vcf_to_pgen.pvar
    }

    meta {
        description: "This workflow is to convert imputation SNP vcf format to pgen format"

    }
}

task convert_vcf_to_pgen {
    input {
        File vcf 
        File vcf_index
        String out_prefix
    }
    command <<<
        plink2 --vcf ~{vcf} --make-pgen --out ~{output_prefix}
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        memory: "6G"
        disks: "local-disk 400 SSD"
    }

    output {
        File pgen = "${out_prefix}.pgen"
        File psam = "${out_prefix}.psam"
        File pvar = "${out_prefix}.pvar"
    }
}