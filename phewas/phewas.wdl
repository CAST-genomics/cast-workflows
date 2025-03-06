version 1.0

workflow phewas {
    input {
        Array[String] chrom
        Array[Int] start
        Array[String] tr_vcf
        Array[String] tr_vcf_index
        File phecode_file
        String GOOGLE_PROJECT
        String WORKSPACE_BUCKET
        Int? mem
        Int? cpu
    }
    
    scatter (i in range(2)) {
        call phetk {
        input:
            chrom=chrom[i],
            start=start[i],
            tr_vcf=tr_vcf[i],
            tr_vcf_index=tr_vcf_index[i],
            phecode_file=phecode_file,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            WORKSPACE_BUCKET=WORKSPACE_BUCKET,
            mem=mem,
            cpu=cpu,
        }
    }

    call merge_outputs {
        input:
            batch_results=phetk.results,
            batch_significant_hits=phetk.significant_hits
    }

    output {
        File phewas_results = merge_outputs.results
        File significant_hits = merge_outputs.significant_hits
    }

    meta {
      description: "Run phewas using phetk on a set of loci."
    }
}

task merge_outputs {
    input {
        Array[File] batch_results
        Array[File] batch_significant_hits
    }

    command <<<
        cat ~{sep=' ' batch_significant_hits} > significant_hits.txt
        cat ~{sep=' ' batch_results} > phewas_results.csv
    >>>

    runtime {
        docker: "ubuntu:latest"
        preemptible: 1
    }

    output {
        File significant_hits="significant_hits.txt"
        File results="phewas_results.csv"
    }
}


task phetk {
    input {
        String chrom
        Int start
        File tr_vcf
        File tr_vcf_index
        File phecode_file
        String GOOGLE_PROJECT
        String WORKSPACE_BUCKET
        Int? mem
        Int? cpu
    }

    String locus_name = "~{chrom}_~{start}"

    command <<<
        export GCS_REQUESTER_PAYS_PROJECT="~{GOOGLE_PROJECT}"
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        export WORKSPACE_BUCKET="~{WORKSPACE_BUCKET}"

        cd /cast-workflows/phewas
        python3 run_phewas.py --locus ~{locus_name} \
                     --chrom ~{chrom} \
                     --start ~{start} \
                     --n-threads ~{cpu} \
                     --no-plot \
                     --tr-vcf ~{tr_vcf} \
                     --outdir . \
                     --phecode-filename ~{phecode_file} \
                     --print-significant-hits > significant_hits.txt
    >>>
    runtime {
        docker:"sarajava/phewas:v1.1"
        memory:  "~{mem} GB"
        cpu: "~{cpu}"
        disks: "local-disk ~{mem} SSD"
        preemptible: 1
    }

    output {
        File results = "~{locus_name}_phewas_results.csv"
        File significant_hits = "significant_hits.txt"
    }

}
