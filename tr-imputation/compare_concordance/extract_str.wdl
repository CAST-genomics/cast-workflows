version 1.0

workflow tr_extraction {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        File str 
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        
    }

    ### Call each chromosome vcf ###
    Int num = length(vcfs)
    scatter (i in range(num)) {
        File vcf = vcfs[i]
        File vcf_index = vcfs_index[i]
        call extract_str{
            input:
                vcf=vcf,
                vcf_index=vcf_index,
                str=str,
                out_prefix=out_prefix,
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN
        }

    }   

    call merge_outputs {
        input:
            vcfs=extract_str.outvcf,
            vcfs_index=extract_str.outvcf_index,
            out_prefix=out_prefix
    } 

    output {
        File outfile = merge_outputs.outvcf
        File outfile_index = merge_outputs.outvcf_index
    }

    meta {
        description: "This workflow extracts selected STRs from AoU imputation vcfs and writes output."
    }
}

task extract_str {
    input {
        File vcf
        File vcf_index
        File str
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        
    }
    String chrom_outprefix = basename(vcf, "annotated.vcf.gz")

    command <<<
        set -e
        echo "the job started"
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        bcftools view -R ~{str} -f"%CHROM\t%POS\t%REF\t%ALT\t[%TGT\t]\n" ~{vcf} -Oz -o "~{out_prefix}_~{chrom_outprefix}.vcf.gz" 
        tabix -p vcf "~{out_prefix}_~{chrom_outprefix}.vcf.gz"
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }
    output {
        File outvcf = "${out_prefix}_${chrom_outprefix}.vcf.gz"
        File outvcf_index = "${out_prefix}_${chrom_outprefix}.vcf.gz.tbi"
    }
}

task merge_outputs {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
    }

    command <<<
        set -e
        bcftools merge --force-samples ~{sep=" " vcfs} -Oz -o ~{out_prefix}_merged.vcf.gz 
        tabix -p vcf ~{out_prefix}_merged.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        preemptible: 1
    }

    output {
        File outvcf = "${out_prefix}_merged.vcf.gz"
        File outvcf_index = "${out_prefix}_merged.vcf.g.tbi"
    }
}