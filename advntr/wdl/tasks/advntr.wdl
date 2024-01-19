version 1.0

workflow run_advntr {

    input {
        String bam_file
        String bam_index
        String output_dir
    }

    call advntr {
        input :
            bam_file = bam_file,
            bam_index = bam_index,
            output_dir = output_dir
    }

    output {
        File log_file = advntr.log_file
        File filtering_out_file = advntr.filtering_out_file
        File keywords_file = advntr.keywords_file
        File unmapped_file = advntr.unmapped_file
        File genotype_output = advntr.genotype_output
    }

    meta {
        description: "This workflow calls adVNTR to genotype VNTRs"
    }
}

task advntr {
    input {
        File bam_file
        File bam_index
        String output_dir
    }

    parameter_meta {
        bam_file: {
          description: "Input bam file",
          localization_optional: true,
          stream: true
        }
        bam_index: {
          description: "Input bam index file",
          localization_optional: true,
          stream: true
        }
    }

    # all output files except for the vcf file are generated in the work_dir.
    String work_dir = "~{output_dir}/work_dir"
    String bam_basename = sub(basename(bam_file), ".bam", "")

    String logging = "~{work_dir}/log_~{bam_basename}.bam.log"
    String filtering_out = "~{work_dir}/filtering_out_~{bam_basename}.unmapped.fasta.txt"
    String keywords = "~{work_dir}/keywords_~{bam_basename}.unmapped.fasta.txt"
    String unmapped = "~{work_dir}/~{bam_basename}.unmapped.fasta"
    String genotype_output = "./~{bam_basename}.genotype.vcf"

    # VNTR_db is placed in the docker file. So the path is within the docker image.
    String vntr_db = "/adVNTR-1.5.0/vntr_db/hg38_VNTRs_by_TRF.db"


    # Set VNTR ids for genotyping based on input.
    # Two options right now: VNTR in the ACAN gene or the list of 52 phenotype associated VNTRs.
    #String vids = "$(cat /adVNTR-1.5.0/vntr_db/phenotype_associated_vntrs_comma.txt)"
    String vids = "290964"

    # Mkdir only needed on the genomequery. On the gcloud there is no need to mkdir.
    #mkdir ~{work_dir}

    command <<<
        export TMPDIR=/tmp
        advntr genotype  \
        --alignment_file ~{bam_file} \
        --models ~{vntr_db}  \
        --working_directory ~{work_dir} \
        -vid ~{vids} \
        --outfmt vcf \
        --disable_logging \
        --pacbio > ~{genotype_output}
    >>>

    runtime {
        docker:"sarajava/advntr-1.5.0:db"
    }

    output {
        File log_file = "~{logging}"
        File filtering_out_file = "~{filtering_out}"
        File keywords_file = "~{keywords}"
        File unmapped_file = "~{unmapped}"
        File genotype_output = "~{genotype_output}"
    }
}
