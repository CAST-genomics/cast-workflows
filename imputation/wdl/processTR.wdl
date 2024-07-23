version 1.0

workflow processTR {
    input {
        File vcf
        File vcf_index
        File ref
        File ref_index
        String out_prefix
        File header_file

    }

    call extract_TR {
        input:
            vcf=vcf,
            vcf_index=vcf_index,      
            out_prefix=out_prefix,
            header_file=header_file
    }

    call prep_beagle_vcf {
        input:
            vcf=extract_TR.outvcf,
            vcf_index=extract_TR.outvcf_index,
            ref=ref,
            ref_index=ref_index,
            out_prefix=out_prefix

    }

    output {
        File outfile = prep_beagle_vcf.outvcf
        File outfile_index = prep_beagle_vcf.outvcf_index
    }

    meta {
        description: "Extract TR from beagle output and process file for TR tools"
    }
}


task extract_TR {
    input {
        File vcf
        File vcf_index
        String out_prefix
        File header_file
    }
    command <<<
        bcftools annotate -h ~{header_file} ~{vcf} -Oz -o ~{out_prefix}.vcf.gz
        bcftools view -i 'ID="."' ~{out_prefix}.vcf.gz -Oz -o ~{out_prefix}_TR.vcf.gz
        tabix -p vcf ~{out_prefix}_TR.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outvcf = "${out_prefix}_TR.vcf.gz"
        File outvcf_index = "${out_prefix}_TR.vcf.gz.tbi"

    }    
}
task prep_beagle_vcf {
    input {
        File vcf
        File vcf_index
        File ref
        File ref_index
        String out_prefix
    }
    command <<<
        bash /usr/bin/trtools_prep_beagle_vcf.sh hipstr ~{ref} ~{vcf} ~{out_prefix}_converted_TR.vcf.gz

    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outvcf = "${out_prefix}_converted_TR.vcf.gz"
        File outvcf_index = "${out_prefix}_converted_TR.vcf.gz.tbi"

    }    

}




