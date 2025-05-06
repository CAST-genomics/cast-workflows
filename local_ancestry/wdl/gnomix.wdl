version 1.0

workflow run_gnomix {
    input {
        File vcf
        File model
        String chrom
        String out_prefix
        File chainfile
        File refpanel
        File refpanel_index
    }

    call liftover_vcf {
        input:
            vcf=vcf,
            out_prefix=out_prefix,
            chainfile=chainfile,
            chrom=chrom
    }

    call beagle {
        input:
            vcf=liftover_vcf.outvcf,
            vcf_index=liftover_vcf.outvcf_index,
            out_prefix=out_prefix,
            refpanel=refpanel,
            refpanel_index=refpanel_index
    }

    call filter_and_split_vcf {
        input:
            vcf=beagle.outvcf,
            vcf_index=beagle.outvcf_index,
            out_prefix=out_prefix
    }

    call gnomix {
        input:
            vcf_array=filter_and_split_vcf.outvcf_array,
            vcf_index_array=filter_and_split_vcf.outvcf_index_array,
            model=model,
            chrom=chrom,
            out_prefix=out_prefix
    }

    call merge_gnomix {
        input:
            gnomix_outputs_msp=gnomix.msp_outfile_array,
            out_prefix=out_prefix
    }

    output {
        File msp_outfile = merge_gnomix.msp_outfile
    }

    meta {
        description: "Run beagle and gnomix on a specified subset of samples on a single chromosome"
    }
}

task liftover_vcf {
    input {
        File vcf
        String out_prefix
        File hg38_ref = "gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File hg19_ref = "gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta"
        File chainfile
        String chrom
    }

    command <<<
        set -e # fail if anything doesn't succeed

        # Liftover to hg19
        echo "Liftover"
        bcftools +liftover --no-version -Oz ~{vcf} -- \
              -f ~{hg19_ref} \
              -s ~{hg38_ref} \
              -c ~{chainfile} \
              --reject ~{out_prefix}.reject.bcf \
              --reject-type b \
              --write-src --drop-tags FORMAT/AD > ~{out_prefix}_hg19_lift.vcf.gz

        # Sort liftover file
        echo "sort"
        df -ah
        ls -ltrh
        bcftools sort -T . -o ~{out_prefix}_hg19.vcf.gz -Oz ~{out_prefix}_hg19_lift.vcf.gz
        tabix -p vcf ~{out_prefix}_hg19.vcf.gz

        # Restrict to target chromosome
        echo "restrict chromosome"
        bcftools view -r ~{chrom} ~{out_prefix}_hg19.vcf.gz -Oz -o ~{out_prefix}.filtered.vcf.gz
        tabix -p vcf ~{out_prefix}.filtered.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs-plugins:latest"
        memory: "25GB"
        disks: "local-disk 60 HDD"
    }

    output {
        File outvcf = "${out_prefix}.filtered.vcf.gz"
        File outvcf_index = "${out_prefix}.filtered.vcf.gz.tbi"
    }    
}

task beagle {
    input {
        File vcf
        File vcf_index
        String out_prefix
        File refpanel
        File refpanel_index
    }

    command <<<
        set -e
        java -Xmx25g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{refpanel} \
            impute=true \
            out=~{out_prefix}_phased
        tabix -p vcf ~{out_prefix}_phased.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/beagle:latest"
        memory: "30GB"
    }

    output {
       File outvcf = "${out_prefix}_phased.vcf.gz"
       File outvcf_index = "${out_prefix}_phased.vcf.gz.tbi"
    }
}

task filter_and_split_vcf {
    input {
        File vcf
        File vcf_index
        String out_prefix
        Int batch_size = 100 # smaller files to gnomix doesn't explode
    }

    command <<<
        set -e

        # Split by batch_size
        bcftools query -l ~{vcf} > sample_list.txt
        mkdir batches
        split -l ~{batch_size} sample_list.txt batches/batch

        # Make group set
        for file in batches/batch*
        do
            batchname=$(basename $file)
            cat $file | awk -v"prefix=$batchname" '{print $1 "\t-\t"prefix}'
        done > sample_groups.txt

        # Run bcftools split and index files
        bcftools plugin split ~{vcf} \
             -G sample_groups.txt -Oz -o .
        for f in batch*.vcf.gz
        do
            tabix -p vcf $f
        done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs-plugins:latest"
        preemptible: 2
    }

    output {
       Array[File] outvcf_array = glob("batch*.vcf.gz")
       Array[File] outvcf_index_array = glob("batch*.vcf.gz.tbi")
    }
}

task gnomix {
    input {
        Array[File] vcf_array
        Array[File] vcf_index_array
        File model
        String chrom
        String out_prefix
        Int total = length(vcf_array)
    }

    command <<<
        set -e
        cd /gnomix
        tar -xzvf ~{model}

        VCFARRAY=(~{sep=" " vcf_array}) # Process one vcf file at a time. Trying to save memory...
        for (( c = 0; c < ~{total}; c++ )) # bash array are 0-indexed ;)
        do
            vcf=${VCFARRAY[$c]}
            echo "Processing $vcf..."
            python3 gnomix.py ${vcf} . ~{chrom} False pretrained_gnomix_models/chr~{chrom}/model_chm_~{chrom}.pkl
            cp query_results.msp /cromwell_root/~{out_prefix}_${c}.msp
        done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/gnomix:latest"
        memory: "25GB"
        preemptible: 2
    }

    output {
       Array[File] msp_outfile_array = glob("*.msp")
    }
}

task merge_gnomix {
    input {
        Array[File] gnomix_outputs_msp
        String out_prefix
        Int total = length(gnomix_outputs_msp)
    }

    command <<<
        set -e
        # First merge msp files
        MSPFILEARRAY=(~{sep=' ' gnomix_outputs_msp})
        head -n 1 ${MSPFILEARRAY[0]} > ~{out_prefix}.msp
        cat ${MSPFILEARRAY[0]} | grep -v "^#Subpopulation" | cut -f 1-6 -d$'\t'> fixedcols.msp
        for (( c = 0; c < ~{total}; c++ ))
        do
            cat ${MSPFILEARRAY[$c]} | grep -v "^#Subpopulation" | cut -f 1-6 -d$'\t' --complement > data_${c}.msp
        done
        paste fixedcols.msp data*.msp >> ~{out_prefix}.msp
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        preemptible: 2
    }

    output {
        File msp_outfile = "~{out_prefix}.msp"
    }
}