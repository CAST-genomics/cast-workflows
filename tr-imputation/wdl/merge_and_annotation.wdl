version 1.0

workflow merge_and_annotation {
    input {
        Array[File] vcf_list = []
        Array[File] idx_list = []
        File ref_vcf 
        File ref_vcf_idx
        File region_file
        Int? region_num = 2
        Int merge_mem
        String out_prefix
        Int anno_mem
    }

    call merge_batch {
        input : 
            TR_vcfs = vcf_list, 
            TR_idxs = idx_list, 
            merge_prefix = out_prefix,
            mem = merge_mem
        
    }
    scatter (i in range(region_num)) {
        call annotaTR {
            input :
                vcf = merge_batch.merged_vcf,
                vcf_index =  merge_batch.merged_vcf_idx,
                ref_vcf = ref_vcf, 
                ref_index = ref_vcf_idx, 
                out_prefix = out_prefix,
                region_file = region_file,
                region_id = i,
                mem = anno_mem 
      }
        
    }
    call concat_pgen {
        input :
            pgen_files = annotaTR.pgen,
            pvar_files = annotaTR.pvar,
            psam_files = annotaTR.psam,
            mem = anno_mem,
            merged_prefix = out_prefix
    }
    call concat_vcf {
        input :
            vcf_files = annotaTR.vcf_files,
            vcf_indexes = annotaTR.vcf_indexes,
            merged_prefix = out_prefix,
            mem = anno_mem
    }
    output {
        File annotated_vcfs  = concat_vcf.concated_vcf
        File annotated_vcfs_indexes  = concat_vcf.concated_vcf_index
        Array[File] pgen_files = concat_pgen.merged_pgens
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
        ulimit -n 800000
        bcftools merge ~{sep=' ' TR_vcfs} -Oz -o ~{merge_prefix}.vcf.gz --threads 2
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

task annotaTR {
    input {
        File vcf
        File vcf_index
        File ref_vcf
        File ref_index
        File region_file
        Int region_id
        String out_prefix

        Int mem
    }
    command <<<
        set -e
        ulimit -n 800000
        input_region=$(awk -v num=~{region_id} 'NR==num+1 {print $1":"$2"-"$3}' ~{region_file})
        annotaTR --vcf ~{vcf} \
            --ref-panel ~{ref_vcf} \
            --out ~{out_prefix}_annotated_~{region_id} \
            --vcftype hipstr \
            --outtype pgen vcf \
            --vcf-outtype z \
            --region ${input_region} \
            --dosages beagleap_norm \
            --ignore-duplicates \
            --match-refpanel-on locid \
            --warn-on-AP-error \
            --update-ref-alt
          tabix -p vcf ~{out_prefix}_annotated_~{region_id}.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-annotatr:yli091230"
        memory: mem+"G"
    }
    
    output {
        File pgen = "${out_prefix}_annotated_${region_id}.pgen"
        File pvar = "${out_prefix}_annotated_${region_id}.pvar"
        File psam = "${out_prefix}_annotated_${region_id}.psam"
        File vcf_files = "${out_prefix}_annotated_${region_id}.vcf.gz"
        File vcf_indexes = "${out_prefix}_annotated_${region_id}.vcf.gz.tbi"

    }
}

task concat_pgen {
    input {
        Array[File] pgen_files
        Array[File] pvar_files
        Array[File] psam_files
        Int mem
        String merged_prefix
    }
    command <<<
        set -e
        ulimit -n 800000
        pgen_array=(~{sep=" " pgen_files})
        pgen_prefix=("${pgen_array[@]}%.*")
        pgen_prefix_sorted=($(printf "%s\n" "${pgen_prefix[@]}" | sort -u -V))
        # sorted_vcf_array=($(printf "%s\n" "${vcf_array[@]}" | sort -V))
        for f in ${pgen_prefix_sorted[@]:1}; do 
          echo ${f} >> merge_list.txt
        done
        cat merge_list.txt
        plink2 --pfile ${pgen_prefix_sorted[0]} \
            --merge-list ./merge_list.txt \
            --make-pgen \
            --out ~{merged_prefix}_annotated 
    >>>
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem+"G"
    }
    output {
        Array[File] merged_pgens = glob("${merged_prefix}_annotated*")
    }
}

task concat_vcf {
    input {
        Array[File] vcf_files
        Array[File] vcf_indexes
        String merged_prefix
        Int mem
    }
    # String vcf_files_str = sep(", ", vcf_files)
    command <<<
				set -e
        ulimit -n 800000
        vcf_array=(~{sep=" " vcf_files})
        sorted_vcf_array=($(printf "%s\n" "${vcf_array[@]}" | sort -V))
        echo "start concat all files..."
        cpu_n=$(nproc)
        bcftools concat --threads ${cpu_n} -Oz -o ~{merged_prefix}_annotated.vcf.gz "${sorted_vcf_array[@]}"
        tabix -p vcf ~{merged_prefix}_annotated.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: mem + "GB"
        disk: "local-disk 150 SSD"
        maxRetries: 1
    }

    output {
        File concated_vcf = "${merged_prefix}_annotated.vcf.gz"
        File concated_vcf_index = "${merged_prefix}_annotated.vcf.gz.tbi"
    }

}












