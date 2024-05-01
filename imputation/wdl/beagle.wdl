version 1.0

workflow beagle {
    input {
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
    }
    #call split {
    #    input:
    #      vcf=vcf, 
    #      vcf_index=vcf_index,
    #      out_prefix=out_prefix,
    #      GOOGLE_PROJECT=GOOGLE_PROJECT,
    #}
    call beagle {
        input : 
          vcf=vcf, 
          vcf_index=vcf_index,
          ref_panel=ref_panel, 
          ref_panel_index=ref_panel_index,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
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
      description: "Run Beagle on a single chromesome with default parameters"
    }
}

task split {
    input {
        File vcf
        File vcf_index
        String out_prefix
        String GOOGLE_PROJECT
    } 
    String int_out_file="~{out_prefix}_partial"

    command <<<
      #export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      #export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "getting a subset of the vcf file based on a region (bcftools)"
      date
      bcftools view ~{vcf} chr21:5030500-5030600 > ~{int_out_file}.vcf
      #echo "getting a subset of the vcf file based on a region (tabix)"
      #date
      #tabix -p vcf ~{vcf} chr21:5030500-5030600 > ~{int_out_file}.vcf
      date
      vcf-sort ~{int_out_file}.vcf | bgzip -c > ~{int_out_file}.sorted.vcf.gz
      tabix -p vcf ~{int_out_file}.sorted.vcf.gz
    >>>
    
  
    runtime {
	docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
       File outfile = "~{int_out_file}.sorted.vcf.gz"
       File outfile_idx = "~{int_out_file}.sorted.vcf.gz.tbi"
    }
}

task beagle {
    input {
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
    } 

    command <<<
      #export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      #export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "Java max heap size"
      java -XX:+PrintFlagsFinal -version | grep HeapSize
      echo "system memory"
      grep MemTotal /proc/meminfo | awk '{print $2}'
      echo "running beagle"
      date
      java -Xmx20g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=20 \
            out=~{out_prefix}
    >>>
    
  
    runtime {
        #docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	docker: "sarajava/beagle:v3"
	memory: "40GB"
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}



task sort_index_beagle {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}
