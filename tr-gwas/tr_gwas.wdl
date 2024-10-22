version 1.0

workflow tr_gwas {
    input {
        Array[File] pgens = [] 
        Array[File] psams = [] 
        Array[File] pvars = [] 
        Array[File] phenotypes = []
        Array[File] cohorts = []
        String GOOGLE_PROJECT = ""
      
    }

    ### Separate workflow for each phenotype ###

    scatter (pheno_file in phenotypes) {      
        call convert_phenotype {
            input:
                pheno=pheno_file,
                GOOGLE_PROJECT=GOOGLE_PROJECT
        }

            scatter (cohort in cohorts) {
                call run_tr_gwas {
                    input:
                        pgens=pgens,
                        psams=psams,
                        pvars=pvars,
                        pheno=convert_phenotype.outfile_pheno,
                        covar=convert_phenotype.outfile_covar,
                        samples=cohort,
                        out_prefix="${pheno_file.basename}_${cohort}_gwas"
                       
                }
            }
    }   
    

    output {
        Array[Array[File]] gwas_outputs = run_tr_gwas.outfile
    }

    meta {
        description: "Use plink2 to run genome-wide TR GWAS. It runs GWAS on each chromosome for each phenotype and merges summary statstics across all chromsomes."
    }
}


task convert_phenotype {
    input {
        File pheno
        String GOOGLE_PROJECT = ""
    }

    command <<<
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        python /usr/bin/convert_phenotype_plink.py --phenotype ~{pheno}
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
    }

    output {ß
       File outfile_pheno = "${pheno}_pheno_plink.txt"
       File outfile_covar = "${pheno}_covar_combined.txt"
    }
}

task run_tr_gwas {
    input {
        Array[File] pgens
        Array[File] psams
        Array[File] pvars
        File pheno
        File covar
        File samples
        String out_prefix
        Int total = length(pgens)
    }

    command <<<
        # Run GWAS on each chrom
        PFILEARRAY=(~{sep=" " pgens})
        gwas_outfiles=""
        for (( c = 0; c < ~{total}; c++ )); # bash array are 0-indexed 
        do
            pfile=${PFILEARRAY[$c]}
            chrom_outprefix=$(basename $pfile .pgen)
            plink2 --pfile ${chrom_outprefix} \
               --pheno ~{pheno} \
               --linear \
               --covar ~{covar} \
               --keep ~{samples} \
               --covar-variance-standardize \
               --out ${chrom_outprefix}
            gwas_outfiles="${gwas_outfiles} ${chrom_outprefix}.glm.linear"
        done

        # Concatenate all results
        cat ${gwas_outfiles} | head -n 1 > ~{out_prefix}.tab
        cat ${gwas_outfiles} | grep -v POS >> ~{out_prefix}.tab
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
    }

    output {
       File outfile = "${out_prefix}.tab"
    }
}