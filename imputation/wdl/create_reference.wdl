version 1.0

workflow create_reference {
    input {
        File vntr_vcf
        File vntr_vcf_index
        String snp_vcf  
        String snp_vcf_index
        String region
        String out_prefix
        String GOOGLE_PROJECT = ""
    }
    call merge_vntr_nearby_snps {
        input:
          vntr_vcf=vntr_vcf, 
          vntr_vcf_index=vntr_vcf_index,
          snp_vcf=snp_vcf, 
          snp_vcf_index=snp_vcf_index,
          region=region,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
    }

    output {
        File snp_outfile = merge_vntr_nearby_snps.snp_outfile 
        File snp_outfile_idx = merge_vntr_nearby_snps.snp_outfile_idx
        File vntr_snp_outfile = merge_vntr_nearby_snps.vntr_snp_outfile
        File vntr_snp_outfile_idx = merge_vntr_nearby_snps.vntr_snp_outfile_idx
    }
    meta {
      description: "Run Beagle on a single chromesome with default parameters"
    }
}

task merge_vntr_nearby_snps {
    input {
        File vntr_vcf
        File vntr_vcf_index
        String snp_vcf
        String snp_vcf_index
        String region
        String out_prefix
        String GOOGLE_PROJECT
    } 
    String int_out_file="~{out_prefix}_partial"

    command <<<
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "getting a subset of the vcf file based on a region (bcftools)"
      bcftools view ~{snp_vcf} ~{region} > ~{int_out_file}.vcf
      date
      bcftools sort -o ~{int_out_file}.sorted.vcf -o ~{int_out_file}.vcf 
      bgzip -c ~{int_out_file}.sorted.vcf > ~{int_out_file}.sorted.vcf.gz
      date
      tabix -p vcf ~{int_out_file}.sorted.vcf.gz
      # ACAN VNTR region, a 10MB window
      bcftools concat --allow-overlaps ~{int_out_file}.sorted.vcf.gz ~{vntr_vcf} > vntr_snp.vcf
      bcftools sort -o vntr_snp.sorted.vcf vntr_snp.vcf
      bgzip -c vntr_snp.sorted.vcf > vntr_snp.sorted.vcf.gz
      tabix -p vcf vntr_snp.sorted.vcf.gz
    >>>
    
  
    runtime {
        #docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
       File snp_outfile = "~{int_out_file}.sorted.vcf.gz"
       File snp_outfile_idx = "~{int_out_file}.sorted.vcf.gz.tbi"
       File vntr_snp_outfile = "vntr_snp.sorted.vcf.gz"
       File vntr_snp_outfile_idx = "vntr_snp.sorted.vcf.gz.tbi"
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
            window=5 \
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
        echo "Number of TRs in the vcf file"
        bcftools view -i 'ID="."' ~{basename}.sorted.vcf.gz | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}
