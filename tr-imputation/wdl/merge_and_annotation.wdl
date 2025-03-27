version 1.0

workflow merge_and_annotation {
    input {
        Array[File] vcf_list = []
        Array[File] idx_list = []
        File ref_vcf 
        File ref_vcf_idx
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
    call annotaTR {
        input :
            vcf = merge_batch.merged_vcf,
            vcf_index =  merge_batch.merged_vcf_idx,
            ref_vcf = ref_vcf, 
            ref_index = ref_vcf_idx, 
            out_prefix = out_prefix, 
            mem = anno_mem 
    }
    output {
        Array[File] annotated_files  = annotaTR.annotated_files
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

task annotaTR {
    input {
        File vcf
        File vcf_index
        File ref_vcf
        File ref_index
        String out_prefix
        Int mem
    }
    command <<<
        set -e
        ulimit -n 800000
        annotaTR --vcf ~{vcf} \
            --ref-panel ~{ref_vcf} \
            --out ~{out_prefix}_annotated \
            --vcftype hipstr \
            --outtype pgen vcf \
            --vcf-outtype z \
            --dosages beagleap_norm \
            --ignore-duplicates \
            --match-refpanel-on locid \
            --warn-on-AP-error \
            --update-ref-alt
          tabix -p vcf ~{out_prefix}_annotated.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-annotatr:latest"
        memory: mem+"G"
    }
    
    output {
        Array[File] annotated_files = glob("*${out_prefix}_annotated*")
    }
}
