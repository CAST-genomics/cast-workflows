import "imputation.wdl" as imputation

workflow imputation_batch {
    input {
        String vcf
        String vcf_index
        String ref_panel
        String ref_panel_bref
        String ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String chrom
        String subset_vcf_path
        Boolean skip_subset_vcf
        Int? mem
        Int? window_size
        Int? overlap
        Array[File] samples_files
	    File regions_file
       }
       scatter (i in range(length(samples_files))) {
               File samples_file = samples_files[i]
               call imputation.imputation as imputation {
                    vcf=vcf,
                    vcf_index=vcf_index,
                    ref_panel=ref_panel,
                    ref_panel_bref=ref_panel_bref,
                    ref_panel_index=ref_panel_index,
                    out_prefix=out_prefix,
                    GOOGLE_PROJECT=GOOGLE_PROJECT,
                    chrom=chrom,
                    subset_vcf_path=subset_vcf_path,
                    skip_subset_vcf=skip_subset_vcf,
                    mem=mem,
                    window_size=window_size,
                    overlap=overlap,
                    samples_file=samples_file,
                    regions_file=regions_file
               }
       }
       call merge_outputs {
            input:
                individual_vcfs = imputation.outfile,
                individual_vcf_indexes = imputation.outfile_index,
                mem=mem
        }

        call sort_index {
            input:
                vcf=merge_outputs.merged_vcfs,
                mem=mem
        }

        output {
            File merged_vcf = sort_index.sorted_vcf
            File merged_vcf_index = sort_index.sorted_vcf_index
        }

        meta {
            description: "This workflow calls Beagle imputation in batches to run in parallel"
        }
}

task sort_index {
    input {
        File vcf
        Int? mem
    }
    String out_prefix = "merged_samples.sorted"
    command <<<
        bcftools sort -Oz ~{vcf} > ~{out_prefix}.vcf.gz && tabix -p vcf ~{out_prefix}.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
	    memory: mem + "GB"
        #bootDiskSizeGb: mem
	    disks: "local-disk " + mem + " SSD"
        maxRetries: 2
    }
    output {
        File sorted_vcf = "~{out_prefix}.vcf.gz"
        File sorted_vcf_index = "~{out_prefix}.vcf.gz.tbi"
    }
}

task merge_outputs {
    input {
        Array[File] individual_vcfs
        Array[File] individual_vcf_indexes
        Int? mem
    }

    String out_prefix = "merged_samples"

    command <<<
        #mergeSTR --vcfs ~{sep=',' individual_vcfs} --out ~{out_prefix}
        # TODO: Work with the -m flag
        bcftools merge -Oz ~{sep=',' individual_vcfs} > ~{out_prefix}.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-5.0.1:latest"
	    memory: mem + "GB"
        #bootDiskSizeGb: mem
	    disks: "local-disk " + mem + " SSD"
        maxRetries: 2
    }
    output {
        File merged_vcfs = "~{out_prefix}.vcf.gz"
    }
}
