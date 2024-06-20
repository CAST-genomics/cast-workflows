version 1.0

workflow imputation {
    input {
        String vcf 
        String vcf_index 
        String ref_panel
        String ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        String chrom
        Int? mem 
        Int? window_size 
        Int? overlap
        File samples_file 
	File regions_file
    }

    call subset_vcf {
    input:
        samples_file=samples_file,
        regions_file=regions_file,
        ref_panel=ref_panel, 
        ref_panel_index=ref_panel_index,
        vcf=vcf,
        vcf_index=vcf_index,
        GOOGLE_PROJECT=GOOGLE_PROJECT,
        GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
        out_prefix=out_prefix
    }
    
    call index_vcf {
        input:
            vcf=subset_vcf.outfile
    }
    call beagle {
        input : 
          vcf=index_vcf.outvcf, 
          vcf_index=index_vcf.outvcf_index,
          ref_panel=ref_panel, 
          ref_panel_index=ref_panel_index,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
          chrom=chrom,
          mem=mem,
          window_size=window_size,
          overlap=overlap
    }
    call sort_index_beagle {
        input :
            vcf=beagle.outfile
    }
    output {
        File outfile = sort_index_beagle.outvcf 
        File outfile_index = sort_index_beagle.outvcf_index
    }
    meta {
      description: "Run Beagle on a subset of samples on a single chromesome with default parameters"
    }
}


task subset_vcf {
    input {
        String vcf
        String vcf_index
        String ref_panel
        String ref_panel_index
        File? samples_file
	File? regions_file
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        String out_prefix=out_prefix
    }

    command <<<
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        # To get a list of sample ids
        #bcftools query -l ~{vcf} > sample_ids.txt
        # The "bcftools head" command was to check the header for the labeling if contigs e.g. chr21 vs 21.
        # bcftools head ~{vcf} > header.txt
        # Subsetting region for each chromesome
        # Select reference samples from the actual samples (to exclude later)
        bcftools query -l ~{ref_panel} > ref_sample_ids.txt
        # Exclude reference samples from the query. Otherwise beagle will give an error.
        grep -v -x -f ref_sample_ids.txt ~{samples_file} > samples_file_clean.txt
        bcftools view -Oz -R ~{regions_file} -S samples_file_clean.txt ~{vcf} > ~{out_prefix}.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: "60GB"
    }

    output {
        File outfile = "${out_prefix}.vcf.gz"
        File ref_sample_ids = "ref_sample_ids.txt"
    }    
}
task index_vcf {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf")

        #bgzip -c ~{vcf}> ~{basename}.vcf.gz && tabix -p vcf ~{basename}.vcf.gz
    command <<<
        tabix -p vcf ~{vcf}
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
        File outvcf = "${basename}.vcf.gz"
        File outvcf_index = "${basename}.vcf.gz.tbi"
    }
}

task beagle {
    input {
        File vcf 
        File vcf_index 
        String ref_panel
        String ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        String chrom
        Int? mem 
        Int? window_size
        Int? overlap
    } 

    command <<<

        #export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        #export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            overlap=~{overlap} \
            chrom=~{chrom} \
            out=~{out_prefix}_output
    >>>
    
    #file upto 300mb use mem=25
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	    memory: mem + "GB"
    }

    output {
       File outfile = "${out_prefix}_output.vcf.gz"
    }
}

task sort_index_beagle {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
        echo "Number of TRs in the genotyped file"
        bcftools view -i 'ID="."' ~{basename}.sorted.vcf.gz | grep -v "^#" | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}
