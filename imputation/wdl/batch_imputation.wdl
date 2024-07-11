version 1.0

import "imputation.wdl" as imputation_t
#import "processTR.wdl" as processTR_t
#import "processSNP.wdl" as processSNP_t
#import "merge_TR_batch.wdl" as merge_TR_batch_t
#import "merge_SNP_batch.wdl" as merge_SNP_batch_t


workflow batch_imputation {
        input {
                String vcf 
                String vcf_index 
                File ref_panel
                String out_prefix
                String GOOGLE_PROJECT = ""
                Int? mem 
                Int? window_size 
                String? region
                Boolean subset_region = false
                Boolean beagle_region = false
                Array[File] samples = []
                Int? disk
                Int? overlap
                File map

        }
    ### Call subsetting samples with batches ###
        Int num_batches = length(samples)
        scatter (i in range(num_batches)) {
                File sample_batch = samples[i]
                call imputation_t.run_imputation as run_imputation {
                    input:
                        sample=sample_batch,
                        vcf=vcf,
                        vcf_index=vcf_index,
                        ref_panel=ref_panel,
                        region=region,
                        GOOGLE_PROJECT=GOOGLE_PROJECT,
                        subset_region=subset_region,
                        beagle_region=beagle_region,
                        out_prefix=out_prefix+".BATCH"+i,
                        mem=mem,
                        window_size=window_size,
                        disk=disk,
                        overlap=overlap,
                        map=map
                }
                ## extract TR from batches of beagle output
                call extract_TR {

                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        out_prefix=out_prefix+".BATCH"+i

                }
        }

        ## use MergeSTR to merge TR
        call  merge_TR_batch {
            input:
                vcfs=extract_TR.outvcf,
                vcfs_index=extract_TR.outvcf_index,
                out_prefix=out_prefix,
                disk=disk
        }

        output {
            File outfile = merge_TR_batch.outfile
            File outfile_index = merge_TR_batch.outfile_index
        }

        meta {
            description: "This workflow run imputation on batches of sample, extract TRs and merge across a single chromosome with default parameters "
        }
                   
}

task extract_TR {
    input {
        File vcf
        File vcf_index
        String out_prefix
    }
    command <<<
        bcftools view -i 'ID="."' ~{vcf} -Oz -o ~{out_prefix}_TR.vcf.gz
        tabix -p vcf ~{out_prefix}_TR.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outvcf = "${out_prefix}_TR.vcf.gz"
        File outvcf_index = "${out_prefix}_TR.vcf.gz.tbi"

    }    
}

task merge_TR_batch {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
        Int? disk
    }
    command <<<
        bcftools merge ~{sep='' vcfs} -Oz -o ${out_prefix}_TR_merged.vcf.gz
        tabix -p vcf ~{out_prefix}_TR_merged.vcf.gz
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        disks: "local-disk ~{disk} SSD"
    }
    output {
        File outfile = "${out_prefix}_TR_merged.vcf.gz"
        File outfile_index = "${out_prefix}_TR_merged.vcf.gz.tbi"
    }

}

