version 1.0

workflow splitVCF {
    input {
        File vcf_file
        String chromosome
        Int sample_size
        Int split_mem
    } 
    call run_split {
        input :
            combined_vcf = vcf_file,
            mem = split_mem,
            sample_size = sample_size,
            outprefix = chromosome
    }

    output {
      Array[File] batched_vcfs = run_split.batched_vcfs
      Array[File] batched_vcf_idxs = run_split.batched_vcf_idxs
      File sample_list = run_split.sample_list
      
    }
}


task run_split {
    input {
        File combined_vcf
        Int mem
        Int sample_size
        String outprefix
    }

    command <<<
    set -e
    ulimit -n 800000
    echo "extract samples first"
    # split by batches
    bcftools query -l ~{combined_vcf} | \
      awk -v group_size=~{sample_size} -v prefix=~{outprefix} 'BEGIN {FS=OFS="\t"} {print $1,"-",prefix"_batch_"int(NR / group_size)+1}' > sample_list.txt
    
    mkdir ./split_by_samples_hg38
    bcftools plugin split ~{combined_vcf} -G sample_list.txt -Oz -o ./split_by_samples_hg38/
    for f in ./split_by_samples_hg38/*.vcf.gz; do tabix -p vcf $f; done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem+"GB"
        disk:"local-disk 200 SSD"
        maxRetries: 1
    }

    output {
      File sample_list = "sample_list.txt"
      Array[File] batched_vcfs = glob("./split_by_samples_hg38/*.vcf.gz")
      Array[File] batched_vcf_idxs = glob("./split_by_samples_hg38/*.vcf.gz.tbi")
    }

}
