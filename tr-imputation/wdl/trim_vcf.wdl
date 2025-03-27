version 1.0

workflow trimpVCF {
    input {
        String outprefix
        Array[File] pvcf_batches
        String qc_thresholds
        String rm_tags 
        Int threads_num
        Int bcftools_mem
    }
    call run_batch {
        input :
            out_prefix = outprefix,
            vcf_list = pvcf_batches,
            qc_inputs = qc_thresholds,
            tags_to_rm = rm_tags,
            n_threads = threads_num,
            mem = bcftools_mem
    }    

    output {
        File final_vcf = run_batch.batch_final_vcf
        File final_vcf_indx = run_batch.batch_final_vcf_idx    
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
    command <<<
        set -e
        ulimit -n 800000
        echo "start job"
        file_idx=1 
        vcf_list_str=(~{sep=" " vcf_list})
        if [[ "~{tags_to_rm}" != "NA" ]] && [[ "~{qc_inputs}" != "NA" ]]
        then
          for pvcf in ${vcf_list_str[@]}; do
            echo "processing id_${file_idx} ${pvcf} filtering both"
            trimmed_vcf=~{out_prefix}_${file_idx}_trimmed.vcf.gz
            bcftools annotate -x "~{tags_to_rm}" ${pvcf} | bcftools norm -m - | bcftools filter -i "~{qc_inputs}" -Oz -o ${trimmed_vcf} --threads ~{n_threads}
            tabix -p vcf ${trimmed_vcf}
            ((file_idx++))
          done
        elif [[ "~{tags_to_rm}" != "NA" ]] && [[ "~{qc_inputs}" == "NA" ]]
        then
          for pvcf in ${vcf_list_str[@]}; do
            echo "processing ${pvcf}"
            trimmed_vcf=~{out_prefix}_${file_idx}_trimmed.vcf.gz
            bcftools annotate -x "~{tags_to_rm}" ${pvcf} | bcftools norm -m - -Oz -o ${trimmed_vcf} --threads ~{n_threads}
            tabix -p vcf ${trimmed_vcf}
            ((file_idx++))
            rm ${trimmed_vcf}
          done
        else
          for pvcf in ${vcf_list_str[@]}; do
            echo "processing ${pvcf}"
            trimmed_vcf=~{out_prefix}_${file_idx}_trimmed.vcf.gz
            bcftools norm -m - ${pvcf} | bcftools filter -i "~{qc_inputs}" -Oz -o ${trimmed_vcf} --threads ~{n_threads}
            tabix -p vcf ${trimmed_vcf}
            ((file_idx++))
            rm ${trimmed_vcf}
          done
        fi

        echo "start combine"
        ls -lh ./*_trimmed.vcf.gz
        file_list=($(ls -d ./*_trimmed.vcf.gz | sort -V))
        echo "-------------"
        bcftools concat --threads ~{n_threads} "${file_list[@]}" -Oz -o ~{out_prefix}_trimmed.vcf.gz  
        tabix -p vcf ~{out_prefix}_trimmed.vcf.gz        
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        disk: "local-disk 150 SSD"
        maxRetries: 1
    }

    output {
        File batch_final_vcf = "${out_prefix}_trimmed.vcf.gz" 
        File batch_final_vcf_idx = "${out_prefix}_trimmed.vcf.gz.tbi"
    }

    meta {
        description: "This is modified on https://github.com/drarwood/vcf_trimmer/tree/master?tab=readme-ov-file"

    }
}
