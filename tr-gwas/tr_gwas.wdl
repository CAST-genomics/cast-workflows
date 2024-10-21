version 1.0

workflow tr_gwas {
    input {
        Array[File] pfamily = [] 
        Array[File] phenotypes = []
        Array[File] samples = []
        String out_prefix    
    }

    ### Looping through each phenotype in the list and convert into PLINK-friendly format###
    scatter (phenotype in phenotypes) {
        File pheno_file = phenotype.rstrip('_phenocovar.csv')
        call convert_phenotype {
            input:
                pheno= pheno_file
        }
        
        ### Looping through each chromosome for each phenotype ###
        scatter (s in range(1,23)) {
            File pgen = "chr"+s+"_annotated.pgen"
            File pvar = "chr"+s+"_annotated.pvar"
            File psam = "chr"+s+"_annotated.psam"

            ###Looping through ALL, EUR and AFR cohort###
            scatter (sample in samples){
                call run_tr_gwas{
                    input:
                        pgen=pgen,
                        pvar=pvar,
                        psam=psam,
                        chrom = "chr"+s+"_annotated",
                        pheno=convert_phenotype.outfile_pheno,
                        covar=convert_phenotype.outfile_covar,
                        sample=sample,
                        out_prefix=pheno_file+"_chr"+s+sample+"_gwas"
                }
            }
            ###concatenating all chromosomes together 
            call concat_gwas_result{
                input:
                    gwas=run_tr_gwas.outfile,
                    out_prefix=pheno_file+"_genome_wide_"+sample+"_gwas"
            }
        }
    }

    output {
        File outfile =  concat_gwas_result.outfile
    }

    meta {
        description: "This workflow is to use plink2 to run genome wide TR GWAS. It run gwas on each chromosome for each phenotype and merge summary statstics across all chromsomes."
    }
}


task convert_phenotype {
    input {
        File pheno
    }

    command <<<
        bash /usr/bin/convert_phenotype_plink.py --phenotype ~{pheno}
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
    }

    output {
       File outfile_pheno = "${pheno}_pheno_plink.txt"
       File outfile_covar = "${pheno}_covar_combined.txt"
    }
}

task run_tr_gwas {
    input {
        File pgen
        File pvar
        File psam                
        File pheno
        File covar
        File sample
        String chrom
        String out_prefix
    }

    command <<<
        plink2 --pfile ~{chrom} \
               --pheno ~{pheno} \
               --linear \
               --covar ~{covar} \
               --keep ~{sample} \
               --covar-variance-standardize \
               --out ~{out_prefix}
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
    }

    output {
       File outfile = "${out_prefix}.phenotype.glm.linear"
    }
}

task concat_gwas_result {
    input {
        Array[File] gwas 
        String out_prefix
    }

    command <<<
        # Extract the header 
        grep "^#" ~{gwas[0]}  > ~{out_prefix}.tab
        # Append rest of gwas result 
        for i in ~{gwas}
        do
            grep -v "^#" $i >> ~{out_prefix}.tab
        done
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/plink2:latest"
        preemptible: 1
    }

    output {
       File outfile = "${out_prefix}.tab"
    }

}