version 1.0

workflow extract_str {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        File str 
        String out_prefix
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
                out_prefix=out_prefix
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
    }
    String chrom_outprefix = basename(vcf, "annotated.vcf.gz")

    command <<<
        set -e
        
        bcftools view -R ~{str} -f"%CHROM\t%POS\t%REF\t%ALT\t[%TGT\t]\n" ~{vcf} -Oz -o "~{out_prefix}_~{chrom_outprefix}.vcf.gz"
        tabix -p vcf "~{out_prefix}_~{chrom_outprefix}.vcf.gz"
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }
    output {
        File outvcf = "${out_prefix}_~{chrom_outprefix}.vcf.gz"
        File outvcf_index = "${out_prefix}_~{chrom_outprefix}.vcf.gz.tbi"
    }
}

task merge_outputs {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
    }

    command <<<
        bcftools merge --force-samples ~{sep=" " vcfs} -Oz -o ~{out_prefix}_merged.vcf.gz
        tabix -p vcf ~{out_prefix}_merged.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outvcf = "~{out_prefix}_merged.vcf.gz"
        File outvcf_index = "~{out_prefix}_merged.vcf.g.tbi"
    }
}