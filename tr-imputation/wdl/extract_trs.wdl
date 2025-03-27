version 1.0

workflow extract_and_merge {
    input {
        Array[File] vcf_list = []
        Array[File] idx_list = []
        Int extract_mem
        String merge_prefix
        Int merge_mem
    }

    Int batch_size = length(vcf_list)
    scatter(i in range(batch_size)) {
        File vcf = vcf_list[i]
        File idx = idx_list[i]
        String curr_file_str = basename(vcf, ".vcf.gz") 

        call extract_variant {
            input :
                vcf = vcf,
                vcf_idx = idx,
                out_prefix = curr_file_str,
                mem = extract_mem
        }
    }

    call merge_batch {
        input : 
            TR_vcfs = extract_variant.imputed_TR,
            TR_idxs = extract_variant.imputed_idx,
            merge_prefix = merge_prefix,
            mem = merge_mem
        
    }

    output {
        File vcf_files = merge_batch.merged_vcf
        File index_files = merge_batch.merged_vcf_idx
    }
}

task extract_variant {
    input {
        File vcf
        File vcf_idx
        String out_prefix
        Int mem
    }
    
    command <<<
        set -e
        ulimit -n 800000

        echo "start extracting TRs..." 
        bcftools view -i 'ID ~"EnsTR"' ~{vcf} -Oz -o ~{out_prefix}_imputed_TR.vcf.gz
        echo "finish extracting, start index now" 
        tabix -p vcf ~{out_prefix}_imputed_TR.vcf.gz
    >>>
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem+"G"
        maxRetries: 1
    }
    output {
        File imputed_TR = "${out_prefix}_imputed_TR.vcf.gz"
        File imputed_idx = "${out_prefix}_imputed_TR.vcf.gz.tbi"
    }
}

task merge_batch {
    input {
        Array[File] TR_vcfs
        Array[File] TR_idxs
        String merge_prefix
        Int mem
    }
    
    command <<<
        bcftools merge ~{sep=' ' TR_vcfs} -Oz -o ~{merge_prefix}.vcf.gz
        tabix -p vcf ~{merge_prefix}.vcf.gz
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem+"G"
        maxRetries: 1
    }
    
    output {
        File merged_vcf = "${merge_prefix}.vcf.gz"
        File merged_vcf_idx = "${merge_prefix}.vcf.gz.tbi"
    }
}
