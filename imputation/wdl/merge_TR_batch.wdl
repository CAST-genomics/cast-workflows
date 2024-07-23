version 1.0

workflow merge_batch {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
        Int? disk
    }

    call mergeSTR {
        input:
            vcfs=vcfs,
            vcfs_index=vcfs_index,
            out_prefix=out_prefix,
            disk=disk
    }

    #call index_vcf {
    #    input:
    #        vcf=mergeSTR.outfile
    #}

    output {
        File outfile = mergeSTR.outfile
        #File outfile_index = index_vcf.outvcf_index
    }

    meta {
        description: "Merge batches of vcf file across a single chromosome"
    }
}

task mergeSTR {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
        Int? disk
        }

    command <<<
        mergeSTR --vcfs ~{sep=',' vcfs} --out ~{out_prefix}_TR_merged --vcftype hipstr
        
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-6.0.1:latest"
        #disks: "local-disk 20 SSD"
        disks: "local-disk ~{disk} SSD"
    }
    output {
        File outfile = "${out_prefix}_TR_merged.vcf"
    }

}

#task index_vcf {
#    input {
#      File vcf
#    }
#    String basename = basename(vcf, ".vcf")
#
#    command <<<
#       vcf-sort ~{vcf} | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
#      
#    >>>
#
#    runtime {
#        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
#    }
#
#    output {
#    File outvcf = "${basename}.sorted.vcf.gz"
#    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
#    }
#}
