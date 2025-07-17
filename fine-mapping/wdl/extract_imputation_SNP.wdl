version 1.0

workflow extract_SNP {
    input {
        Array[File] vcf = [] 
        Array[File] vcf_index = []  
        String out_prefix 
    }

    call extract_SNP {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            out_prefix=out_prefix
    }

    call merge_batch {
        input:
            vcfs=vcfs,
            vcfs_index=vcfs_index,
            out_prefix=out_prefix

    }
    output {
        File outfile_vcf = merge_batch.outvcf
        File outfile_vcf_index = merge_batch.outvcf_index
    }

    meta {
        description: "This workflow is to extract SNP from Imputation file per chromosome and merge batches "

    }
}

task extract_SNP {
    input {
        File vcf = [] 
        File vcf_index = []  
        String out_prefix 
    }
    command <<<
        set -e
        bcftools view -i 'ID!~"EnsTR"' ~{vcf} -Oz -o ~{out_prefix}_SNP.vcf.gz
        tabix -p vcf ~{out_prefix}_SNP.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        preemptible: 1
    }

    output {
        File outvcf = "${out_prefix}_SNP.vcf.gz"
        File outvcf_index = "${out_prefix}_SNP.vcf.gz.tbi"

    }
}

task merge_batch {
    input {
        Array[File] vcfs = []
        Array[File] vcfs_index = []
        String out_prefix 
    }

    command <<<
        bcftools merge ~{sep=' ' vcfs} -Oz -o ~{out_prefix}_SNP_merged.vcf.gz 
        tabix -p vcf ~{out_prefix}_SNP_merged.vcf.gz        
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: "25G"
        disks: "local-disk 300 SSD"
    }
    output {
        File outvcf = "${out_prefix}_SNP_merged.vcf.gz"
        File outvcf_index = "${out_prefix}_SNP_merged.vcf.gz.tbi"
    }

}

