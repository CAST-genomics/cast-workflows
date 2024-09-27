version 1.0

workflow beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        String out_prefix
        String GOOGLE_PROJECT = ""
        Int? mem 
        Int? window_size
        Boolean beagle_region = false
        String? region
        Int? overlap
        File map
        Int? disk
    }


    call beagle {
        input : 
          vcf=vcf,
          vcf_index=vcf_index,
          ref_panel=ref_panel, 
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          mem=mem,
          window_size=window_size,
          beagle_region=beagle_region,
          region=region,
          map=map,
          overlap=overlap,
          disk=disk
    }

    output {
        File outfile = beagle.outvcf 
        File outfile_index = beagle.outvcf_index
    }
    meta {
      description: "Run Beagle on a subset of samples on a single chromesome with default parameters"
    }
}

task beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        String out_prefix
        String GOOGLE_PROJECT = ""
        Int? mem 
        Int? window_size
        String? region
        Int? overlap
        File map
        Int? disk
        Boolean beagle_region = false
   
    } 

    command <<<

        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
         
        if [[ "~{beagle_region}" == false ]] ; then
        java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            overlap=~{overlap} \
            out=~{out_prefix}_output \
            map=~{map} 
        else
            java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            chrom=~{region} \
            overlap=~{overlap} \
            out=~{out_prefix}_output \
            map=~{map} 
        fi

        tabix -p vcf ~{out_prefix}_output.vcf.gz
        #gsutil cp ~{out_prefix}_output.vcf.gz ~{out_prefix}_output.vcf.gz.tbi $WORKSPACE_BUCKET/imputation_result/~{out_prefix}
        
    >>>
    
    #file upto 300mb use mem=25
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	    memory: mem + "GB"
        disks: "local-disk ${disk} SSD"
    }

    output {
       File outvcf = "${out_prefix}_output.vcf.gz"
       File outvcf_index = "${out_prefix}_output.vcf.gz.tbi"
    }
}