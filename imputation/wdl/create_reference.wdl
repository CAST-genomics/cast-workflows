version 1.0

workflow create_reference {
    input {
        String vntr_vcf
        String vntr_vcf_index
        String snp_vcf  
        String snp_vcf_index
        String region
        String samples
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
          samples=samples,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
    }

    output {
        File snp_outfile = merge_vntr_nearby_snps.snp_outfile 
        File snp_outfile_idx = merge_vntr_nearby_snps.snp_outfile_idx
        File vntr_snp_vcf = merge_vntr_nearby_snps.vntr_snp_vcf
        File vntr_snp_outfile = merge_vntr_nearby_snps.vntr_snp_outfile
        File vntr_snp_np_outfile = merge_vntr_nearby_snps.vntr_snp_np_outfile
        File vntr_snp_outfile_idx = merge_vntr_nearby_snps.vntr_snp_outfile_idx
    }
    meta {
      description: "Run Beagle on a single chromesome with default parameters"
    }
}

task merge_vntr_nearby_snps {
    input {
        String vntr_vcf
        String vntr_vcf_index
        String snp_vcf
        String snp_vcf_index
        String region
        String samples
        String out_prefix
        String GOOGLE_PROJECT
    } 
    String int_out_file="~{out_prefix}_partial"
    String vntr_vcf_sample_subset="vntr_vcf_sample_subset.vcf.gz"
    #String samples_file="gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/saraj/vntr_reference_panel/vntr_samples_clean_sorted.txt"

    command <<<
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "getting a subset of the vcf file based on a region and a sample list(bcftools)"
      bcftools view -O z -r ~{region} -S ~{samples} -e 'POS=87296520' ~{snp_vcf} > ~{int_out_file}.vcf.gz
      date
      bcftools sort -O z  ~{int_out_file}.vcf.gz > ~{int_out_file}.sorted.vcf.gz
      ls -lht
      echo "number of variants in unsorted file"
      zcat ~{int_out_file}.vcf.gz | grep -v "^#" | wc -l
      echo "number of variants in sorted file"
      zcat ~{int_out_file}.sorted.vcf.gz | grep -v "^#" | wc -l
      rm ~{int_out_file}.vcf.gz
      date
      tabix -p vcf ~{int_out_file}.sorted.vcf.gz
      #bcftools view -h ~{int_out_file}.sorted.vcf.gz > header.txt
      # ACAN VNTR region, a 10MB window
      echo "subset samples from the vntr file"
      bcftools view -O z -S ~{samples} ~{vntr_vcf} > ~{vntr_vcf_sample_subset}
      tabix -p vcf ~{vntr_vcf_sample_subset}
      echo "concat with sorted files"
      bcftools concat -O z --allow-overlaps ~{int_out_file}.sorted.vcf.gz ~{vntr_vcf_sample_subset} > vntr_snp.vcf.gz
      echo "sort the resulting file"
      bcftools sort -O z vntr_snp.vcf.gz > vntr_snp.sorted.vcf.gz
      echo "number of variants in pre sorted file"
      zcat vntr_snp.vcf.gz | grep -v "^#" | wc -l
      zcat vntr_snp.sorted.vcf | grep -v "^#" | wc -l
      #bgzip -c vntr_snp.sorted.vcf > vntr_snp.sorted.vcf.gz
      tabix -p vcf vntr_snp.sorted.vcf.gz

      #echo "concat with unsorted files"
      echo "empty" > vntr_snp_no_pre_sorted.vcf.gz
      #bcftools concat --allow-overlaps ~{int_out_file}.vcf ~{vntr_vcf} > vntr_snp_no_pre_sorted.vcf
      #echo "number of variants in no pre sorted file"
      #cat vntr_snp_no_pre_sorted.vcf | grep -v "^#" | wc -l
      #bgzip -c vntr_snp_no_pre_sorted.vcf > vntr_snp_no_pre_sorted.vcf.gz
      #bcftools sort -o vntr_snp.sorted.vcf vntr_snp.vcf
      ls -lth
    >>>
    
  
    runtime {
        memory: "60GB"
        #docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
       File snp_outfile = "~{int_out_file}.sorted.vcf.gz"
       File snp_outfile_idx = "~{int_out_file}.sorted.vcf.gz.tbi"
       File vntr_snp_vcf = "vntr_snp.vcf.gz"
       File vntr_snp_np_outfile = "vntr_snp_no_pre_sorted.vcf.gz"
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
