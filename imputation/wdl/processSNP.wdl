version 1.0

workflow processSNP {
    input {
        File vcf
        File vcf_index
        String out_prefix
    }

    call extract_SNP {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            out_prefix=out_prefix
    }


    output {
        File outfile = extract_SNP.outvcf
        File outfile_index = extract_SNP.outvcf_index
    }

    meta {
        description: "Extract and merge SNP from beagle output"
    }
}

task extract_SNP {
    input {
        File vcf
        File vcf_index
        String out_prefix
    }

    command <<<
        bcftools view -i 'ID!="."' ~{vcf} -Oz -o ~{out_prefix}_SNP.vcf.gz
        tabix -p vcf ~{out_prefix}_SNP.vcf.gz

    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outvcf = "${out_prefix}_SNP.vcf.gz"
        File outvcf_index = "${out_prefix}_SNP.vcf.gz.tbi"

    }    
}

