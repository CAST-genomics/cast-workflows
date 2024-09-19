version 1.0

import "create_reference.wdl" as create_ref

workflow create_reference_batch {
    input {
        String vntr_vcf
        String vntr_vcf_index
        String snp_vcf  
        String snp_vcf_index
        String regions
        String samples
        String out_prefix
        String GOOGLE_PROJECT = ""
        Int window
        Int mem
        String map
    }

    scatter (chrom_idx in [21, 22]) {
        String chrom="chr~{chrom_idx}"
        String map_file="~{map}/plink.~{chrom}_w_chr.GRCh38.map"
        call create_ref.create_reference as create_reference {
          input:
            vntr_vcf=vntr_vcf,
            vntr_vcf_index=vntr_vcf_index,
            snp_vcf=snp_vcf,
            snp_vcf_index=snp_vcf_index,
            regions=regions,
            samples=samples,
            out_prefix=out_prefix,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            window=window,
            mem=mem,
            map=map_file,
            chrom=chrom,
        }
    }

    call concat {
      input:
        vcfs=create_reference.phased_vntr_snp_vcf,
        vcf_indexes=create_reference.phased_vntr_snp_index,
        mem=mem,
    }

    call sort_index {
      input:
        vcf=concat.merged_vcfs,
	mem=mem*2
    }

    call create_ref.bref as bref {
      input:
        vcf=sort_index.sorted_vcf,
        vcf_index=sort_index.sorted_vcf_index,
        mem=mem,
    }

    output {
        File merged_ref_vcf = sort_index.sorted_vcf
        File merged_ref_vcf_index = sort_index.sorted_vcf_index
        File merged_ref_bref = bref.outfile   
    }

    meta {
      description: "Create reference for the whole genome"
    }

}

task sort_index {
    input {
        File vcf
        Int? mem
    }
    String out_prefix = "merged_ref.sorted"
    command <<<
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

task concat {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        Int? mem
    }

    String out_prefix = "merged_ref_files"

    command <<<
        bcftools concat -Oz ~{sep=' ' vcfs} > ~{out_prefix}.vcf.gz
    >>>
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
