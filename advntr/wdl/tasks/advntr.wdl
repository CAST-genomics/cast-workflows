version 1.0

workflow run_advntr {

    input {
        File bam_file
        String bam_file_advntr_path
        File bam_index
        String output_dir
        String sample_id
    }

    call advntr {
        input :
            bam_file = bam_file,
            bam_file_advntr_path = bam_file_advntr_path,
            bam_index = bam_index,
            output_dir = output_dir,
            sample_id = sample_id,
    }

    output {
        File? log_file = advntr.log_file
        File? filtering_out_file = advntr.filtering_out_file
        File? keywords_file = advntr.keywords_file
        File? unmapped_file = advntr.unmapped_file
        File genotype_output = advntr.genotype_output
    }

    meta {
        description: "This workflow calls adVNTR to genotype VNTRs"
    }
}


task advntr {
    input {
        File bam_file
        String bam_file_advntr_path
        File bam_index
        String output_dir
        String sample_id
    }

    parameter_meta {
        bam_file: {
          description: "Input bam file",
        }
        bam_index: {
          description: "Input bam index file",
        }
    }

    # all output files except for the vcf file are generated in the work_dir.
    String work_dir = "."

    String logging = "~{work_dir}/log_~{sample_id}.bam.log"
    String filtering_out = "~{work_dir}/filtering_out_~{sample_id}.unmapped.fasta.txt"
    String keywords = "~{work_dir}/keywords_~{sample_id}.unmapped.fasta.txt"
    String unmapped = "~{work_dir}/~{sample_id}.unmapped.fasta"
    String genotype_output = "./~{sample_id}.genotype.vcf"

    # VNTR_db is placed in the docker file. So the path is within the docker image.
    String vntr_db = "/adVNTR/vntr_db/hg38_VNTRs_by_TRF.db"


    # Set VNTR ids for genotyping based on input.
    # Two options right now: VNTR in the ACAN gene or the list of 52 phenotype associated VNTRs.
    #String vids = "$(cat /adVNTR/vntr_db/phenotype_associated_vntrs_comma.txt)"
    String vids = "290964"

    command <<<
        export TMPDIR=/tmp
        ls -lh ~{work_dir}
        ls -lh ~{work_dir}/fc-secure*/*/*
        /usr/bin/time -v advntr genotype  \
        --alignment_file ~{bam_file_advntr_path} \
        --models ~{vntr_db}  \
        --working_directory ~{work_dir} \
        -vid ~{vids} \
        --outfmt vcf \
        --pacbio > ~{genotype_output}
        ls -lh ~{work_dir}
    >>>

    runtime {
        docker:"sarajava/advntr:1.5.0_db_u18"
        cpu: "4"
        memory: "18G"
    }

    output {
        File? log_file = "~{logging}"
        File? filtering_out_file = "~{filtering_out}"
        File? keywords_file = "~{keywords}"
        File? unmapped_file = "~{unmapped}"
        File genotype_output = "~{genotype_output}"
    }
}
