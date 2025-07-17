version 1.0

workflow extract_SNP {
    input {
        Array[File] vcfs = [] 
        Array[File] vcfs_index = []  
        String out_prefix 
    }


    Int num_batches = length(vcfs)
    scatter (i in range(num_batches)) {
        File batch = vcfs[i]
        File batch_index = vcfs[i]+".tbi"
            call extract_SNP {
                input:
                    vcfs=batch,
                    vcfs_index=batch_index,
                    out_prefix=out_prefix+".BATCH"+i
            }
    }
    call merge_batch {
        input:
            vcfs=extract_SNP.outvcf,
            vcfs_index=extract_SNP.outvcf_index,
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
        File vcfs 
        File vcfs_index 
        String out_prefix 
    }
    command <<<
        set -e
        bcftools view -i 'ID!~"EnsTR"' ~{vcfs} -Oz -o ~{out_prefix}_SNP.vcf.gz
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
        Array[File] vcfs 
        Array[File] vcfs_index 
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

