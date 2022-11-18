version 1.0

workflow run_advntr {

    input {
        File bam = "/nucleus/projects/saraj/vntr/data/advntr_simulated_test_data/CSTB_2_5_testdata.bam"
        File bam_index = "/nucleus/projects/saraj/vntr/data/advntr_simulated_test_data/CSTB_2_5_testdata.bam.bai"
        File vntr_db = "/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/hg38_selected_VNTRs_Illumina.db"
    }

    call advntr {
        input :
            bam = bam,
            bam_index = bam_index,
            vntr_db = vntr_db,
    }

    output {
        File log_file = advntr.log_file
        File filtering_out_file = advntr.filtering_out_file
        File keywords_file = advntr.keywords_file
        File unmapped_file = advntr.unmapped_file
    }

    meta {
        description: "This workflow calls adVNTR to genotype VNTRs"
    }
}

task advntr {
    input {
        File bam
        File bam_index
        File vntr_db
    }

    String working_directory = "./working_directory"
    String bam_basename = sub(basename(bam), ".bam", "")

    String logging = "~{working_directory}/log_~{bam_basename}.bam.log"
    String filtering_out = "~{working_directory}/filtering_out_~{bam_basename}.unmapped.fasta.txt"
    String keywords = "~{working_directory}/keywords_~{bam_basename}.unmapped.fasta.txt"
    String unmapped = "~{working_directory}/~{bam_basename}.unmapped.fasta"

    command <<<
        echo "Working_directory: ~{working_directory}"
        echo "Input bam file: ~{bam}"
        echo "VNTR database: ~{vntr_db}"
        echo "Output Log file: ~{logging}"
        echo "Output filtering file: ~{filtering_out}"
        echo "Output keywords file: ~{keywords}"
        echo "Output unmapped file: ~{unmapped}"
        mkdir ~{working_directory}
        advntr genotype \
            --alignment_file ~{bam} \
            --models ~{vntr_db} \
            --working_directory ~{working_directory} \
            --disable_logging
    >>>

    runtime {
        docker:"sarajava/advntr:1.5.0"
    }

    output {
        File log_file = "~{logging}"
        File filtering_out_file = "~{filtering_out}"
        File keywords_file = "~{keywords}"
        File unmapped_file = "~{unmapped}"
    }
}
