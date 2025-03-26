version 1.0

workflow bfileTovcf {
    input {
        File bed
        File bim
        File fam
        File ref_genome
        Int convert_mem 
        String chrom
    }

    call convert {
        input :
            bed = bed,
            bim = bim,
            fam = fam,
            ref_genome = ref_genome,
            outprefix = chrom,
            mem = convert_mem
    }

    output {
        File converted_vcf = convert.combined_vcf
    }
}

task convert {

    input {
        File bed
        File bim
        File fam
        File ref_genome
        String outprefix
        Int mem
    }

    command <<<
        set -e
        echo "start conversion"
        
        mkdir -p bgenToVCF
        plink2 --bed ~{bed} \
            --bim ~{bim} \
            --fam ~{fam} \
            --export vcf id-paste=iid bgz \
            --fa ~{ref_genome} \
            --ref-from-fa \
            --out "bgenToVCF/~{outprefix}"
        
        echo "Adding information"
        tabix -p vcf bgenToVCF/~{outprefix}.vcf.gz

        bcftools +fill-tags "bgenToVCF/~{outprefix}.vcf.gz" -Oz -o "bgenToVCF/~{outprefix}_extra_tags_added_hg19.vcf.gz" -- -t AF,AN,AC
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem + "GB" 
        maxRetries: 1 
    }

    output {
        File combined_vcf = "bgenToVCF/~{outprefix}_extra_tags_added_hg19.vcf.gz"
    }
}
