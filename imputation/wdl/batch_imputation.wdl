version 1.0

import "beagle.wdl" as imputation_t
import "processTR.wdl" as processTR_t
#import "processSNP.wdl" as processSNP_t
import "merge_TR_batch.wdl" as merge_TR_batch_t
#import "merge_SNP_batch.wdl" as merge_SNP_batch_t


workflow batch_imputation {
        input {
                File ref_vcf 
                File ref_index 
                File ref_panel
                String out_prefix
                String GOOGLE_PROJECT = ""
                Int? mem 
                Int? merge_mem
                Int? window_size 
                String? region
                Boolean beagle_region = false
                Array[File] sample_batches = []
                Int? disk
                Int? overlap
                File map
                

        }
    ### Call subsetting samples with batches ###
        
        Int num_batches = length(sample_batches)
        scatter (i in range(num_batches)) {
                File batch = sample_batches[i]
                File batch_index = sample_batches[i]+".tbi"
                call imputation_t.beagle as run_imputation {
                    input:
                        vcf=batch,
                        vcf_index= batch_index,
                        ref_panel=ref_panel,
                        region=region,
                        GOOGLE_PROJECT=GOOGLE_PROJECT,
                        beagle_region=beagle_region,
                        out_prefix=out_prefix+".BATCH"+i,
                        mem=mem,
                        window_size=window_size,
                        disk=disk,
                        overlap=overlap,
                        map=map
                }
                ## extract TR from batches of beagle output
                call processTR_t.processTR as processTR {
                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        out_prefix=out_prefix+".BATCH"+i
     
                }
                ## extract SNP from batches of beagle output
                #call processSNP_t.processSNP as processSNP {
                #    input:
                #        vcf=run_imputation.outfile,
                #        vcf_index=run_imputation.outfile_index,
                #        out_prefix=out_prefix+".BATCH"+i

                #}
        }

        ## use bcftools to merge TR
        call merge_TR_batch_t.merge_batch as merge_TR_batch {
            input:
                vcfs=processTR.outfile,
                vcfs_index=processTR.outfile_index,
                out_prefix=out_prefix,
                disk=disk,
                merge_mem=merge_mem
        }

        ## use bcftools to merge SNP
        #call merge_SNP_batch_t.merge_SNP_batch as merge_SNP_batch {
        #    input:
        #        vcfs=processSNP.outfile,
        #        vcfs_index=processSNP.outfile_index,
        #        out_prefix=out_prefix

        #}

        call annotaTR {
            input:
                vcf=merge_TR_batch.outfile,
                vcf_index=merge_TR_batch.outfile_index,
                ref_vcf=ref_vcf,
                ref_index=ref_index,
                out_prefix=out_prefix,
                merge_mem=merge_mem,
                disk=disk

        }

        output {
            File outfile_pgen = annotaTR.pgen
            File outfile_psam = annotaTR.psam
            File outfile_pvar = annotaTR.pvar
            
        }

        meta {
            description: "This workflow run imputation on batches of sample, extract TRs, merge  across a single chromosome and run annotaTR with default parameters "
        }
                   
}


task annotaTR {
    input {
        File vcf 
        File vcf_index
        File ref_vcf
        File ref_index
        String out_prefix
        Int? merge_mem
        Int? disk
    }
    
    command <<<
        annotaTR --vcf ~{vcf} \
                 --ref-panel ~{ref_vcf} \
                 --out ~{out_prefix}_annotated \
                 --vcftype hipstr \
                 --outtype pgen \
                 --dosages bestguess_norm \
                 --ignore-duplicates \
                 --match-refpanel-on locid \
                 --fix-bcftools-offset
        

    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-6.0.2:latest"
        disks: "local-disk ~{disk} SSD"
        memory: merge_mem + "GB"
    }

    output {
        File pgen = "${out_prefix}_annotated.pgen"
        File psam = "${out_prefix}_annotated.psam"
        File pvar = "${out_prefix}_annotated.pvar"
    }
}