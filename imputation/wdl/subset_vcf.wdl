version 1.0

workflow subset_vcf {
    input {
        String vcf 
        String vcf_index 
        String out_prefix
        File samples_file 
        File regions_file
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        
    }

    call subset_vcf {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            out_prefix=out_prefix,
            samples_file=samples_file,
            regions_file=regions_file,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN
            
    }

    call index_vcf {
        input:
            vcf=subset_vcf.outfile
    }

    output {
        File outfile = index_vcf.outvcf 
        File outfile_index = index_vcf.outvcf_index
    }
    meta {
      description: "Subset samples on a single chromesome with default parameters"
    }
}

task subset_vcf {
    input {
        String vcf
        String vcf_index
        String out_prefix=out_prefix
        File samples_file
        File regions_file
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""    
    }


    command <<<
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        bcftools view -S ~{samples_file} ~{vcf} > ~{out_prefix}.vcf

    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
       File outfile = "${out_prefix}.vcf"
    }
}

task index_vcf {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf")

    command <<<
        bgzip -c ~{vcf}> ~{basename}.vcf.gz && tabix -p vcf ~{basename}.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
        File outvcf = "${basename}.vcf.gz"
        File outvcf_index = "${basename}.vcf.gz.tbi"
    }
}
