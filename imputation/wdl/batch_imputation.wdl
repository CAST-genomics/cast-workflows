version 1.0

import "imputation.wdl" as imputation_t
import "processTR.wdl" as processTR_t
#import "processSNP.wdl" as processSNP_t
import "merge_TR_batch.wdl" as merge_TR_batch_t
#import "merge_SNP_batch.wdl" as merge_SNP_batch_t


workflow batch_imputation {
        input {
                String vcf 
                String vcf_index 
                File ref_panel
                File ref
                File ref_index
                String out_prefix
                String GOOGLE_PROJECT = ""
                Int? mem 
                Int? window_size 
                String? region
                Boolean subset_region = false
                Boolean beagle_region = false
                Array[File] samples = []
                File header_file
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
                call processTR_t.processTR as processTR {
                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        ref=ref,
                        ref_index=ref_index,
                        out_prefix=out_prefix+".BATCH"+i,
                        header_file=header_file

                }
                ## extract SNP from batches of beagle output
                #call processSNP_t.processSNP as processSNP {
                #    input:
                #        vcf=run_imputation.outfile,
                #        vcf_index=run_imputation.outfile_index,
                #        out_prefix=out_prefix+".BATCH"+i

                #}
        }

        ## use MergeSTR to merge TR
        call merge_TR_batch_t.merge_batch as merge_TR_batch {
            input:
                vcfs=processTR.outfile,
                vcfs_index=processTR.outfile_index,
                out_prefix=out_prefix,
                disk=disk
        }

        ## use bcftools to merge SNP
        #call merge_SNP_batch_t.merge_SNP_batch as merge_SNP_batch {
        #    input:
        #        vcfs=processSNP.outfile,
        #        vcfs_index=processSNP.outfile_index,
        #        out_prefix=out_prefix

        #}


        output {
            File trvcf = merge_TR_batch.outfile
            Array[File] outvcf = run_imputation.outfile
            Array[File] outvcf_index = run_imputation.outfile_index
            #File trvcf_index = merge_TR_batch.outfile_index
            #File snpvcf = merge_SNP_batch.outfile
            #File snpvcf_index = merge_SNP_batch.outfile_index
        }

        meta {
            description: "This workflow run imputation on batches of sample, extract TRs and merge across a single chromosome with default parameters "
    }
                   
}

