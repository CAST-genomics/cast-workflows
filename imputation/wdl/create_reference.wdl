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
   
    call subset_snps {
        input:
          snp_vcf=snp_vcf, 
          snp_vcf_index=snp_vcf_index,
          all_regions=regions,
          samples=samples,
          out_prefix=out_prefix,
          chrom=chrom,
          mem=mem,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
    }
    call merge_vntr_snps {
        input:
          snp_vcf=subset_snps.snp_vcf, 
          snp_vcf_index=subset_snps.snp_vcf_index,
          vntr_vcf=vntr_vcf, 
          vntr_vcf_index=vntr_vcf_index,
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
          chrom=chrom,
          vcf=beagle_phase.outfile,
          mem=mem,
    }
    call add_tags {
        input:
          vcf=sort_index_beagle.outvcf,
          vcf_index=sort_index_beagle.outvcf_index,
          annotation_vcf=vntr_vcf,
          annotation_vcf_index=vntr_vcf_index,
          mem=mem,
    }
    call bref {
        input:
          mem=mem,
          vcf=add_tags.outvcf,
          vcf_index=add_tags.outvcf_index,
          chrom=chrom,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
    }
    output {
        File phased_vntr_snp_vcf = add_tags.outvcf
        File phased_vntr_snp_index = add_tags.outvcf_index
        File phased_vntr_snp_bref = bref.outfile
    }
    meta {
      description: "Create reference based on genotyped VNTRs and nearby SNPs"
    }
}

task subset_snps {
    input {
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
      df -h
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "Copying regions file"
      gsutil cp ~{all_regions} .
      grep -w ~{chrom} ~{all_regions_basename} > ~{regions}
      echo "Region file (current chromosome)"
      cat ~{regions}
      echo "getting a subset of the vcf file based on a region and a sample list(bcftools)"
      bcftools view -O z -I -R ~{regions} ~{snp_vcf} > ~{snp_region_file}_no_sample_sort.vcf.gz
      echo "Sorting sample ids"
      bcftools view -O z -S ~{samples} ~{snp_region_file}_no_sample_sort.vcf.gz > ~{snp_region_file}.vcf.gz
      bcftools sort -O z -o ~{snp_region_file}.sorted.vcf.gz ~{snp_region_file}.vcf.gz
      echo "number of variants in unsorted file"
      zcat ~{snp_region_file}.vcf.gz | grep -v "^#" | wc -l
      date
      tabix -p vcf ~{snp_region_file}.vcf.gz
    >>>
    
  
    runtime {
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        preemptible: 1
    }

    output {
       File snp_vcf = "~{snp_region_file}.vcf.gz"
       File snp_vcf_index = "~{snp_region_file}.vcf.gz.tbi"
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
    String regions="regions.bed"
    String all_regions_basename=basename(all_regions)
    String vntr_vcf_sample_subset="vntr_vcf_sample_subset.vcf.gz"

    command <<<
      df -h
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      echo "Copying regions file"
      gsutil cp ~{all_regions} .
      grep -w ~{chrom} ~{all_regions_basename} > ~{regions}
      echo "Region file (current chromosome)"
      cat ~{regions}
      echo "subset samples from the vntr file"
      bcftools view -O z -S ~{samples} -R ~{regions} ~{vntr_vcf} > ~{vntr_vcf_sample_subset}
      tabix -p vcf ~{vntr_vcf_sample_subset}
      echo "concat with sorted files"
      bcftools concat -O z --allow-overlaps ~{snp_vcf} ~{vntr_vcf_sample_subset} > vntr_snp.vcf.gz
      echo "sort the resulting file"
      bcftools sort -O z -o vntr_snp.sorted.vcf.gz vntr_snp.vcf.gz
      echo "number of variants in pre sorted file"
      zcat vntr_snp.vcf.gz | grep -v "^#" | wc -l
      echo "number of variants in sorted file"
      zcat vntr_snp.sorted.vcf | grep -v "^#" | wc -l
      echo "number of VNTRs in sorted file"
      zcat vntr_snp.sorted.vcf | grep -v "^#" | awk '$3 ~ /chr[0-9]*_[0-9]*$/{print}' | wc -l
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
        String GOOGLE_PROJECT
        String chrom
    }
    String basename = basename(vcf, ".vcf.gz")
    command <<<
        zcat ~{vcf} | java -jar /bref3.jar > ~{basename}.bref3
        gsutil -u ~{GOOGLE_PROJECT} cp ~{vcf} gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/saraj/vntr_reference_panel/p_g_vntrs/phased/~{chrom}/
        gsutil -u ~{GOOGLE_PROJECT} cp ~{vcf_index} gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/saraj/vntr_reference_panel/p_g_vntrs/phased/~{chrom}/
        gsutil -u ~{GOOGLE_PROJECT} cp ~{basename}.bref3 gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/saraj/vntr_reference_panel/p_g_vntrs/phased/~{chrom}/
    >>>
    runtime {
        docker:"sarajava/beagle:v4"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
        preemptible: 1
    }
    output {
        File outfile="~{basename}.bref3"
    }
}

task sort_index_beagle {
    input {
      File vcf
      String chrom
      Int mem
    }

    String outfile_tmp="vntr_ref_~{chrom}.sorted.tmp.vcf.gz"
    String outfile="vntr_ref_~{chrom}.sorted.vcf.gz"

    command <<<
        set -e
        tabix -p vcf ~{vcf}
        bcftools sort -Oz ~{vcf} > ~{outfile_tmp} && tabix -p vcf ~{outfile_tmp}
        bcftools view -h ~{outfile_tmp} | head -n 4 > header.txt
        echo "##source=adVNTR ver. 1.5.0" >> header.txt
        bcftools view -h ~{outfile_tmp} | tail -n +5 >> header.txt
        bcftools reheader -h header.txt ~{outfile_tmp} | bcftools view -Oz > ~{outfile}
        tabix -p vcf ~{outfile}
        echo "Number of TRs in the vcf file"
        bcftools view  ~{outfile} | grep -v "^#" | awk '$3 ~ /chr[0-9]*_[0-9]*$/{print}' | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
	memory: mem + "GB"
	disks: "local-disk " + mem + " SSD"
    }

    output {
    File outvcf = "~{outfile}"
    File outvcf_index = "~{outfile}.tbi"
  }
}

task add_tags {
    input {
      File vcf
      File vcf_index
      File annotation_vcf
      File annotation_vcf_index
      Int mem
    }

   String basename = basename(vcf, ".vcf.gz")
   String outfile="~{basename}.annotate.vcf.gz"

   command <<<
       bcftools annotate -Oz -a ~{annotation_vcf} -c CHROM,POS,VID,RU ~{vcf} > ~{outfile}
       tabix -p vcf ~{outfile}
   >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
	memory: mem + "GB"
	disks: "local-disk " + mem*4 + " SSD"
        preemptible: 1
    }

    output {
    File outvcf = "~{outfile}"
    File outvcf_index = "~{outfile}.tbi"
  }
}
