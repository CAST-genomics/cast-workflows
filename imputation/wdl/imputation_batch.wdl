version 1.0

import "imputation.wdl" as imputation_workflow

workflow imputation_batch {
    input {
        String vcf
        String vcf_index
        String ref_panel
        String ref_panel_bref
        String ref_panel_index
        String genetic_map
        String out_prefix
        String GOOGLE_PROJECT
        String chrom
        #String subset_vcf_path
        Boolean skip_subset_vcf
        Int vid
        String motif
        Int? mem
        Int? window_size
        Int? overlap
        Array[File] samples_files
        File regions_file
       }
       scatter (i in range(246)) {
               #File samples_file = samples_files[i]
               String vcf_file="~{vcf}/chr15_batch~{i+1}.vcf.gz"
               String vcf_index_file="~{vcf_file}.tbi"
               File samples_file=samples_files[0]
               call imputation_workflow.imputation as imputation {
                 input:
                    vcf=vcf_file,
                    vcf_index=vcf_index_file,
                    ref_panel=ref_panel,
                    ref_panel_bref=ref_panel_bref,
                    ref_panel_index=ref_panel_index,
                    genetic_map=genetic_map,
                    out_prefix=out_prefix,
                    GOOGLE_PROJECT=GOOGLE_PROJECT,
                    chrom=chrom,
                    #subset_vcf_path=subset_vcf_path,
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
                vid=vid,
                motif=motif,
                mem=mem*2
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
        Int vid
        String motif
        Int? mem
    }
    String out_prefix = "merged_samples.sorted"
    String sorted_prefix = "merged_samples_pre_reheader"
    command <<<
        echo "Sorting vcf file"
        bcftools sort -Oz ~{vcf} > ~{sorted_prefix}.vcf.gz && tabix -p vcf ~{sorted_prefix}.vcf.gz
        echo "Updating the header"
        bcftools view -h ~{sorted_prefix}.vcf.gz | grep "^##" > header.txt
        echo "##source=adVNTR ver. 1.5.0" >> header.txt
        echo '##INFO=<ID=VID,Number=1,Type=Integer,Description="VNTR id in the VNTR database">' >> header.txt
        echo '##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit or consensus motif of the VNTR">' >> header.txt
        bcftools view -h ~{sorted_prefix}.vcf.gz | grep -v "^##" >> header.txt
        bcftools reheader -h header.txt ~{sorted_prefix}.vcf.gz > ~{sorted_prefix}_rh.vcf.gz
        echo "Adding VID info field"
        zcat ~{sorted_prefix}_rh.vcf.gz | sed 's/END=/VID=~{vid};RU=~{motif};END=/g' > ~{sorted_prefix}_w_vid.vcf
        bcftools view -Oz ~{sorted_prefix}_w_vid.vcf > ~{out_prefix}.vcf.gz
        tabix -p vcf ~{out_prefix}.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
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
        touch ~{sep=' ' individual_vcfs}
        bcftools merge -Oz ~{sep=' ' individual_vcfs} > ~{out_prefix}.vcf.gz
    >>>
        # TODO: Work with the -m flag
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
        maxRetries: 2
    }
    output {
        File merged_vcfs = "~{out_prefix}.vcf.gz"
    }
}
