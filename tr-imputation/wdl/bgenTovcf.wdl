version 1.0

workflow bgenTovcf {
    input {
        File bgen_file
        File sample_file
#        File ref_genome
        Int batch_size
        Int split_mem 
        String chrom
        String bgen_qc
    }

    call convert_split {
        input :
            input_file = bgen_file,
            input_sample = sample_file,
#            ref_genome = ref_genome,
            outprefix = chrom,
            sample_size = batch_size, 
            mem = split_mem,
            qc_option = bgen_qc
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
        File input_file
        File input_sample
#        File ref_genome
        String outprefix
        Int sample_size
        Int mem
        String qc_option
    }

    command <<<
        set -e
        echo "start conversion"
        
        mkdir -p bgenToVCF
        plink2 --bgen ~{input_file} "ref-first"\
            --sample ~{input_sample} ~{qc_option} \
            --export vcf id-paste=iid bgz \
            --out "bgenToVCF/~{outprefix}"
        
        echo "Adding information"
        tabix -p vcf bgenToVCF/~{outprefix}.vcf.gz

        bcftools +fill-tags "bgenToVCF/~{outprefix}.vcf.gz" -Oz -o "bgenToVCF/~{outprefix}_extra_tags_added_hg19.vcf.gz" -- -t AF,AN,AC
        ls -lht ./bgenToVCF 
        echo "finish converting, start splitting..." 
        rm ./bgenToVCF/~{outprefix}.vcf.gz
        
        echo "extract samples first" 
        # split by batches
        bcftools query -l ./bgenToVCF/~{outprefix}_extra_tags_added_hg19.vcf.gz | \
            awk -v group_size=~{sample_size} 'BEGIN {FS=OFS="\t"} {print $1,"batch"int(NR / group_size)+1}' > sample_list.txt
        
        echo "extract samples first"
        mkdir batch_files/ 
        awk -v outdir="batch_files/" 'BEGIN {FS=OFS="\t"} {print $1 > outdir$2".txt"}' sample_list.txt
       
        echo "start splitting into $(ls batch_files/*.txt | wc -l) batches"
        mkdir ./split_by_samples_hg19
        idx=1 
        for batch in batch_files/*.txt; do
            batch_name=$(basename ${batch} | cut -d "." -f 1)
            echo "processing ${idx} ${batch}"
            bcftools view -S ${batch} ./bgenToVCF/~{outprefix}_extra_tags_added_hg19.vcf.gz -Oz -o ./split_by_samples_hg19/~{outprefix}_${batch_name}_hg19.vcf.gz
            tabix -p vcf ./split_by_samples_hg19/~{outprefix}_${batch_name}_hg19.vcf.gz
            ((idx++))
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
