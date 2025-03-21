version 1.0

workflow bfileTovcf {
    input {
        File bed
        File bim
        File fam
        File ref_genome
        Int batch_size
        Int split_mem 
        String chrom
    }

    call convert_split {
        input :
            bed = bed,
            bim = bim,
            fam = fam,
            ref_genome = ref_genome,
            outprefix = chrom,
            sample_size = batch_size, 
            mem = split_mem
    }

    output {
        File sample_list = convert_split.sample_list
        File converted_vcf = convert_split.combined_vcf
        Array[File] batched_vcf = convert_split.batched_vcfs 
        Array[File] batched_idxs = convert_split.batched_vcf_idxs
    }
}

task convert_split {

    input {
        File bed
        File bim
        File fam
        File ref_genome
        String outprefix
        Int sample_size
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
        ls -lht ./bgenToVCF 
        echo "finish converting, start splitting..." 
        rm ./bgenToVCF/~{outprefix}.vcf.gz
        
        echo "extract samples first" 
        # split by batches
        bcftools query -l ~{outprefix}_extra_tags_added_hg19.vcf.gz | \
            awk -v group_size=~{sample_size} 'BEGIN {FS=OFS="\t"} {print $1,"batch"int(NR / group_size)+1}' > sample_list.txt
        
        echo "extract samples first"
        mkdir batch_files/ 
        awk -v outdir="batch_files/" 'BEGIN {FS=OFS="\t"} {print $1 > outdir$2".txt"}' sample_list.txt
       
        echo "start splitting into $(ls batch_files/*.txt | wc -l) batches"
        mkdir ./split_by_samples_hg19
        for batch in batch_files/*.txt; do
            batch_name=$(basename ${batch} | cut -d "." -f 1)
            bcftools view -S ${batch} ~{outprefix}_extra_tags_added_hg19.vcf.gz -Oz -o ./split_by_samples_hg19/~{outprefix}_${batch_name}_hg19.vcf.gz
            tabix -p vcf ./split_by_samples_hg19/~{outprefix}_${batch_name}_hg19.vcf.gz
        done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem + "GB" 
        maxRetries: 1 
    }

    output {
        File combined_vcf = "bgenToVCF/~{outprefix}_extra_tags_added_hg19.vcf.gz"
        Array[File] batched_vcfs = glob("./split_by_samples_hg19/*.gz")
        Array[File] batched_vcf_idxs = glob("./split_by_samples_hg19/*.tbi")
        File sample_list = "sample_list.txt"
    }
}

#task liftover {
#      input {
#          Array[File] batched_hg19_vcfs
#      }
#
#      command <<<
#          set -e
#          echo "start liftover"
#          mkdir ./processed_folder
#          hg19_vcf_list=(~{sep=" " batched_hg19_vcfs}) 
#          for hg19_vcf in ${hg19_vcf_list[@]}; do 
#            echo "start processing ${hg19_vcf}" 
#            liftover.py ${hg19_vcf} "./processed_folder"  
#            ls ./processed_folder/${hg19_vcf%.vcf.gz}_hg38.vcf 
#            bgzip ./processed_folder/${hg19_vcf%.vcf.gz}_hg38.vcf
#            # sort ? 
#            rm  ./processed_folder/${hg19_vcf%.vcf.gz}_hg38.vcf
#            bgzip ./processed_folder/${hg19_vcf%.vcf.gz}_unlifted.vcf 
#      >>>
#      
#      runtime {
#          docker: "gcr.io/ucsd-medicine-cast/pyliftover:latest"
#          memory: mem + "GB"
#          disk: 
#      }
#
#      output {
#          Array[File] batched_hg38_vcfs=glob("./processed_folder/*hg38.vcf.gz")
#          Array[File] batched_hg38_idxs=glob("./processed_folder/*hg38.vcf.gz.tbi")
#          Array[File] batched_hg38_unlifted=glob("../processed_folder/*unlifted.vcf.gz")
#          Array[File] liftover_logs=glob("./processed_folder/*liftOver.log")
#      }
#
#}
