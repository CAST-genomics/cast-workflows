version 1.0

workflow beagle {
    input {
	File map
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
    }
    call beagle {
        input : 
	  map=map,
          vcf=vcf, 
          vcf_index=vcf_index,
          ref_panel=ref_panel, 
          ref_panel_index=ref_panel_index,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
    }

    output {
        File outfile = beagle.out_reference
    }
    meta {
      description: "Run Beagle on long read genotype data to generate an imputation reference for long VNTRs"
    }
}


task beagle {
    input {
	File map
        File vcf
        File vcf_index
        File ref_panel
        File ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
    } 

    command <<<
      echo "Java max heap size"
      java -XX:+PrintFlagsFinal -version | grep HeapSize
      echo "system memory"
      grep MemTotal /proc/meminfo | awk '{print $2}'
      echo "running beagle"
      date
      java -Xss5m -Xmx20g \
      	-jar /beagle.jar \
        gt=~{vcf} \
        ref=~{ref_panel} \
        map=~{map} \
        out=~{out_prefix} \
        nthreads=16 \
        niterations=10 \
        ne=20000 \
        impute=true \
        gprobs=true \
        window=20 \
        seed=-99999 

     date
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
