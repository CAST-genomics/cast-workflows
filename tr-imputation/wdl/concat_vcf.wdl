version 1.0

workflow run_concat{
    input {
        Array[File] vcf_files
        Array[File] vcf_indexs
        String merged_prefix
        Int mem
    }

    call concat_vcf {
        input:
            vcf_files = vcf_files,
            vcf_indexs = vcf_indexs,
            merged_prefix = merged_prefix,
            mem = mem
    }

    output {
        File final_vcf = concat_vcf.concated_vcf
#        File final_vcf_index = concat_vcf.concated_vcf_index
    }

}

task concat_vcf {
    input {
        Array[File] vcf_files
        Array[File] vcf_indexs
        String merged_prefix
        Int mem
    }
    # String vcf_files_str = sep(", ", vcf_files) 
    command <<<
        ulimit -n 800000
        vcf_array=(~{sep=" " vcf_files})
        sorted_vcf_array=($(printf "%s\n" "${vcf_array[@]}" | sort -V))
        echo "start concat all files..."
        cpu_n=$(nproc)
        bcftools concat --threads ${cpu_n} -Oz -o ~{merged_prefix}_sorted.vcf.gz "${sorted_vcf_array[@]}"  
        # tabix -p vcf ~{merged_prefix}_sorted.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        disk: "local-disk 150 SSD"
        maxRetries: 1 
    }

    output {
        File concated_vcf = "${merged_prefix}_sorted.vcf.gz"
#        File concated_vcf_index = "${merged_prefix}_sorted.vcf.gz.tbi"
    }

}

