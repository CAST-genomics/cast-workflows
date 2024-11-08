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
        String bucket
        String chrom
        Int? mem
        Int? window_size
        Int? overlap
        Array[File] samples_files
        File regions_file
       }
       scatter (i in range(246)) {
               # Samples file is no longer used because we are taking all the samples
               # in each downloaded batch.
               File samples_file=samples_files[0]
               String vcf_file="~{vcf}/~{chrom}_batch~{i+1}.vcf.gz"
               String vcf_index_file="~{vcf_file}.tbi"
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
                chrom=chrom,
                mem=mem*2
        }

        call add_tags {
            input:
                vcf=sort_index.sorted_vcf,
                vcf_index=sort_index.sorted_vcf_index,
                annotation_vcf=ref_panel,
                annotation_vcf_index=ref_panel_index,
                bucket=bucket,
                chrom=chrom,
                mem=mem,
        }
    
        call annotaTR {
            input:
	            vcf=add_tags.outvcf,
	            vcf_index=add_tags.outvcf_index,
	            ref_vcf=ref_panel,
	            ref_index=ref_panel_index,
                chrom=chrom,
                mem=mem,
        }

        call copy_to_bucket {
            input:
                pgen=annotaTR.pgen,
                psam=annotaTR.psam,
                pvar=annotaTR.pvar,
                bucket=bucket,
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                chrom=chrom,
                mem=mem,
        }

        output {
            File outfile_pgen = copy_to_bucket.pgen_out
            File outfile_psam = copy_to_bucket.psam_out
            File outfile_pvar = copy_to_bucket.pvar_out
        }

        meta {
            description: "This workflow calls Beagle imputation in batches to run in parallel"
        }
}

task sort_index {
    input {
        File vcf
        String chrom
        Int? mem
    }
    String out_prefix = "imputed_~{chrom}.sorted"
    command <<<
        set -e
        echo "Sorting vcf file"
        bcftools sort -Oz ~{vcf} > ~{out_prefix}.vcf.gz && tabix -p vcf ~{out_prefix}.vcf.gz
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
        set -e
        touch ~{sep=' ' individual_vcf_indexes}
        bcftools merge --merge id -Oz ~{sep=' ' individual_vcfs} > ~{out_prefix}.vcf.gz
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	    memory: mem*4 + "GB"
	    disks: "local-disk " + mem*4 + " SSD"
    }
    output {
        File merged_vcfs = "~{out_prefix}.vcf.gz"
    }
}

task copy_to_bucket {
    input {
        File pgen
        File psam
        File pvar
        String bucket
        String GOOGLE_PROJECT
        String chrom
        Int? mem
    }

   String basename = basename(pgen, ".pgen")
    command <<<
        set -e
        ls -lt
        echo ~{pgen}
        gsutil cp ~{pgen} ~{bucket}/saraj/imputation_output/~{chrom}/
        gsutil cp ~{psam} ~{bucket}/saraj/imputation_output/~{chrom}/
        gsutil cp ~{pvar} ~{bucket}/saraj/imputation_output/~{chrom}/
    >>>
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	    memory: mem + "GB"
        preemptible: 1
    }
    output {
        File pgen_out = "~{pgen}"
        File psam_out = "~{psam}"
        File pvar_out = "~{pvar}"
    }
}

task add_tags {
    input {
        File vcf
        File vcf_index
        File annotation_vcf
        File annotation_vcf_index
        String bucket
        String chrom
        Int? mem
    }

   String basename = basename(vcf, ".vcf.gz")
   String outfile="~{basename}.annotated.vcf.gz"

   command <<<
       set -e
       touch ~{vcf_index} ~{annotation_vcf_index}
       bcftools annotate -Oz -a ~{annotation_vcf} -c CHROM,POS,VID,RU ~{vcf} > ~{outfile}
       tabix -p vcf ~{outfile}
       gsutil cp *.vcf.gz* ~{bucket}/saraj/imputation_output/~{chrom}/
   >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	    memory: mem + "GB"
        preemptible: 1
    }

    output {
        File outvcf = "~{outfile}"
        File outvcf_index = "~{outfile}.tbi"
  }
}

task annotaTR {
    input {
        File vcf 
        File vcf_index
        File ref_vcf
        File ref_index
        String chrom
        Int? mem
    }
    
    String basename = basename(vcf, ".vcf.gz")
    command <<<
        set -e
        echo "basename: ~{basename}"
        annotaTR --vcf ~{vcf} \
                 --ref-panel ~{ref_vcf} \
                 --out ~{basename} \
                 --vcftype advntr \
                 --outtype pgen \
                 --dosages bestguess_norm
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-6.0.2:latest"
        disks: "local-disk ~{mem} SSD"
        memory: mem + "GB"
    }

    output {
        File pgen = "~{basename}.pgen"
        File psam = "~{basename}.psam"
        File pvar = "~{basename}.pvar"
    }
}
