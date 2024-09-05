version 1.0

workflow create_reference {
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
        String chrom
    }
    call merge_vntr_snps {
        input:
          vntr_vcf=vntr_vcf, 
          vntr_vcf_index=vntr_vcf_index,
          snp_vcf=snp_vcf, 
          snp_vcf_index=snp_vcf_index,
          all_regions=regions,
          samples=samples,
          out_prefix=out_prefix,
          chrom=chrom,
          mem=mem,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
    }
    call beagle_phase {
        input:
          vcf=merge_vntr_snps.vntr_snp_vcf,
          vcf_index=merge_vntr_snps.vntr_snp_vcf_index,
          mem=mem,
          window=window,
          map=map,
          chrom=chrom,
    }
    call sort_index_beagle {
        input:
          vcf=beagle_phase.outfile
    }
    call bref {
        input:
          mem=mem,
          vcf=sort_index_beagle.outvcf,
          vcf_index=sort_index_beagle.outvcf_index
    }
    output {
        File phased_vntr_snp_vcf = sort_index_beagle.outvcf
        File phased_vntr_snp_index = sort_index_beagle.outvcf_index
        File phased_vntr_snp_bref = bref.outfile
    }
    meta {
      description: "Run Beagle on a single chromesome with default parameters"
    }
}

task merge_vntr_snps {
    input {
        String vntr_vcf
        String vntr_vcf_index
        String snp_vcf
        String snp_vcf_index
        String all_regions
        String samples
        String chrom
        String out_prefix
        String GOOGLE_PROJECT
        Int mem
    } 
    String snp_region_file="~{out_prefix}_partial"
    String regions="regions.bed"
    String all_regions_basename=basename(all_regions)
    String vntr_vcf_sample_subset="vntr_vcf_sample_subset.vcf.gz"

    command <<<
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "Copying regions file"
      gsutil cp ~{all_regions} .
      grep -w ~{chrom} ~{all_regions_basename} > ~{regions}
      echo "Region file (current chromosome)"
      cat ~{regions}
      echo "getting a subset of the vcf file based on a region and a sample list(bcftools)"
      bcftools view -O z -R ~{regions} -S ~{samples} -e 'POS=87296520' ~{snp_vcf} > ~{snp_region_file}.vcf.gz
      echo "Sorting"
      bcftools sort -O z ~{snp_region_file}.vcf.gz > ~{snp_region_file}.sorted.vcf.gz
      echo "number of variants in unsorted file"
      zcat ~{snp_region_file}.vcf.gz | grep -v "^#" | wc -l
      echo "number of variants in sorted file"
      zcat ~{snp_region_file}.sorted.vcf.gz | grep -v "^#" | wc -l
      rm ~{snp_region_file}.vcf.gz
      date
      tabix -p vcf ~{snp_region_file}.sorted.vcf.gz
      # VNTR region
      echo "subset samples from the vntr file"
      bcftools view -O z -S ~{samples} -R ~{regions} ~{vntr_vcf} > ~{vntr_vcf_sample_subset}
      tabix -p vcf ~{vntr_vcf_sample_subset}
      echo "concat with sorted files"
      bcftools concat -O z --allow-overlaps ~{snp_region_file}.sorted.vcf.gz ~{vntr_vcf_sample_subset} > vntr_snp.vcf.gz
      echo "sort the resulting file"
      bcftools sort -O z vntr_snp.vcf.gz > vntr_snp.sorted.vcf.gz
      echo "number of variants in pre sorted file"
      zcat vntr_snp.vcf.gz | grep -v "^#" | wc -l
      zcat vntr_snp.sorted.vcf | grep -v "^#" | wc -l
      tabix -p vcf vntr_snp.sorted.vcf.gz
    >>>
    
  
    runtime {
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
       File vntr_snp_vcf = "vntr_snp.sorted.vcf.gz"
       File vntr_snp_vcf_index = "vntr_snp.sorted.vcf.gz.tbi"
    }
}

task beagle_phase {
    input {
        File vcf
        File vcf_index
        Int mem
        Int window
        File map
        String chrom
    } 
    String basename = basename(vcf, ".vcf.gz")
    String outfile = "phased_~{basename}"

    command <<<
      echo "Java max heap size"
      java -XX:+PrintFlagsFinal -version | grep HeapSize
      echo "system memory"
      grep MemTotal /proc/meminfo | awk '{print $2}'
      echo "running beagle for phasing"
      date
      java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            chrom=~{chrom} \
            map=~{map} \
            window=~{window} \
            out=~{outfile}
    >>>
    
  
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
    }

    output {
       File outfile = "~{outfile}.vcf.gz"
    }
}

task bref {
    input {
        File vcf
        File vcf_index
        Int mem
    }
    String basename = basename(vcf, ".vcf.gz")
    command <<<
        zcat ~{vcf} | java -jar /bref3.jar > ~{basename}.bref3
    >>>
    runtime {
        docker:"sarajava/beagle:v4"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
    }
    output {
        File outfile="~{basename}.bref3"
    }
}

task sort_index_beagle {
    input {
      File vcf
    }

    command <<<
        tabix -p vcf ~{vcf}
        bcftools sort -Oz ~{vcf} > vntr_ref.sorted.vcf.gz && tabix -p vcf vntr_ref.sorted.vcf.gz
        echo "Number of TRs in the vcf file"
        bcftools view -i 'ID="."' vntr_ref.sorted.vcf.gz | grep -v "^#" | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "vntr_ref.sorted.vcf.gz"
    File outvcf_index = "vntr_ref.sorted.vcf.gz.tbi"
  }
}
