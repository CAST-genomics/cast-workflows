version 1.0

workflow trimpVCF {
    input {
        String outprefix
        Array[Array[File]] pvcf_batches=[]
        String qc_thresholds
        String rm_tags 
        Int threads_num=1
        Int bcftools_mem=16
    }
    
    Int num_batches = length(pvcf_batches)
    scatter(i in range(num_batches)) {
        Array[File] pvcfs=pvcf_batches[i]
        call run_batch {
            input : 
                out_prefix = outprefix+"_BATCH"+i,
                vcf_list = pvcfs,
                qc_inputs = qc_thresholds,
                tags_to_rm = rm_tags,
                n_threads = threads_num,
                mem = bcftools_mem
      }
  
  }

    call concat_vcf {
        input :
            vcf_files = run_batch.batch_final_vcf,
            vcf_indexs = run_batch.batch_final_vcf_idx,
            merged_prefix = outprefix + "_ALL",
            mem = bcftools_mem
    }

    output {
        File final_vcf = concat_vcf.concated_vcf
        File final_vcf_indx = concat_vcf.concated_vcf_index
    }
}

task run_batch {
    input {
        String out_prefix
        Array[File] vcf_list 
        String qc_inputs
        String tags_to_rm
        Int n_threads
        Int mem
    }
    # String vcf_list_str = sep(" ", vcf_list)
    command <<<
        echo "start job"
        # cd ~/
        # pwd
        file_idx=1 
        vcf_list_str=(~{sep=" " vcf_list})
        if [[ "~{tags_to_rm}" != "NA" ]] && [[ "~{qc_inputs}" != "NA" ]]
        then
          for pvcf in ${vcf_list_str[@]}; do
            echo "processing id_${file_idx} ${pvcf} filtering both"
            trimmed_vcf=~{out_prefix}_${file_idx}_trimmed.vcf.gz
            bcftools annotate -x "~{tags_to_rm}" ${pvcf} | bcftools norm -m - | bcftools filter -i "~{qc_inputs}" -Oz -o ${trimmed_vcf} --threads ~{n_threads}
            
            sorted_vcf=~{out_prefix}_${file_idx}_trimmed_sorted.vcf.gz
            bcftools sort -Oz -o ${sorted_vcf} ${trimmed_vcf}
            tabix -p vcf ${sorted_vcf}
            rm ${trimmed_vcf} 
            ((file_idx++))
          done
        elif [[ "~{tags_to_rm}" != "NA" ]] && [[ "~{qc_inputs}" == "NA" ]]
        then
          for pvcf in ${vcf_list_str[@]}; do
            echo "processing ${pvcf}"
            trimmed_vcf=~{out_prefix}_${file_idx}_trimmed.vcf.gz
            bcftools annotate -x "~{tags_to_rm}" ${pvcf} | bcftools norm -m - -Oz -o ${trimmed_vcf} --threads ~{n_threads}

            sorted_vcf=~{out_prefix}_${file_idx}_trimmed_sorted.vcf.gz
            bcftools sort -Oz -o ${sorted_vcf} ${trimmed_vcf}
            tabix -p vcf ${sorted_vcf}
            ((file_idx++))
            rm ${trimmed_vcf}
          done
        else
          for pvcf in ${vcf_list_str[@]}; do
            echo "processing ${pvcf}"
            trimmed_vcf=~{out_prefix}_${file_idx}_trimmed.vcf.gz
            bcftools norm -m - ${pvcf} | bcftools filter -i "~{qc_inputs}" -Oz -o ${trimmed_vcf} --threads ~{n_threads}

            sorted_vcf=~{out_prefix}_${file_idx}_trimmed_sorted.vcf.gz
            bcftools sort -Oz -o ${sorted_vcf} ${trimmed_vcf}
            tabix -p vcf ${sorted_vcf}
            ((file_idx++))
            rm ${trimmed_vcf}
          done
        fi
        echo "start combine"
        ls ./*_sorted.vcf.gz
        echo "-------------\n"
        bcftools concat -a  ./*_sorted.vcf.gz -Oz -o ~{out_prefix}.vcf.gz
        echo "sort and index"
        bcftools sort -Oz -o ~{out_prefix}_sorted.vcf.gz ~{out_prefix}.vcf.gz
        tabix -p vcf ~{out_prefix}_sorted.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        cpu: 8
        disks: "local-disk 200 SSD"
        maxRetries: 3
    }

    output {
        File batch_final_vcf = "${out_prefix}_sorted.vcf.gz" 
        File batch_final_vcf_idx = "${out_prefix}_sorted.vcf.gz.tbi"
    }

    meta {
        description: "This is modified on https://github.com/drarwood/vcf_trimmer/tree/master?tab=readme-ov-file"

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
        echo "start concat all files..."
        # vcf_files_str=~{sep=" " vcf_files}
        bcftools concat -a -Oz -o ~{merged_prefix}.vcf.gz ~{sep=" " vcf_files}
        ls 
        bcftools sort -Oz -o ~{merged_prefix}_sorted.vcf.gz ~{merged_prefix}.vcf.gz
        tabix -p vcf ~{merged_prefix}_sorted.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        cpu: 8
        disk: "local-disk 200 SSD"
        maxRetries: 3
    }

    output {
        File concated_vcf = "${merged_prefix}_sorted.vcf.gz"
        File concated_vcf_index = "${merged_prefix}_sorted.vcf.gz.tbi"
    }

}

