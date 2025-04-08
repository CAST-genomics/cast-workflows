version 1.0

workflow splitVCF {
    input {
        File vcf_file
        String chromosome
        Int batch_num
        Int sample_size
        Int split_mem
    }
    call get_batch {
        input :
            combined_vcf = vcf_file,
            outprefix = chromosome,
            sample_size = sample_size,
            number_of_batch = batch_num,
            mem = split_mem
    }
    scatter (batch_file in get_batch.batched_split) {
        call run_split {
            input :
                combined_vcf = vcf_file,
                batch_file = batch_file,
                mem = split_mem
        }

    }

    output {
      Array[Array[File]] batched_vcfs = run_split.batched_vcfs
      Array[Array[File]] batched_vcf_idxs = run_split.batched_vcf_idxs
      Array[File] sample_list = get_batch.batched_split
      
    }
    meta {
        workflow_id: "workflow-GzpbYx0JX3JXF743Jj3pbfgz"
    }
}

task get_batch {
    input {
        File combined_vcf
        String outprefix
        Int sample_size
        Int number_of_batch
        Int mem
    }

    command <<<
        set -e 
        ulimit -n 800000
        echo "split sample files to batches"    
        n_batch=~{number_of_batch}
        bcftools query -l ~{combined_vcf} | \
          awk -v group_size=~{sample_size} -v prefix=~{outprefix} -v batch_n=${n_batch} '
          BEGIN {
              FS=OFS="\t";
              file_idx = 1;
              line_count = 0;
          }
          {
              lines[line_count++] = $0;  # Store all input lines
          }
          END {
              total_batches = int((line_count + group_size - 1) / group_size);  # Compute total batches
              batches_per_file = int((total_batches + (batch_n-1)) / batch_n);  # Distribute batches across 8 files

              for (i = 0; i < line_count; i++) {
                  batch = int(i / group_size) + 1;  # Compute batch number
                  file_idx = int((batch - 1) / batches_per_file) + 1;  # Assign batch to one of 8 files
                  file = prefix"_sample_list_" file_idx ".txt";
                  print lines[i], "-", prefix"_batch_" batch > file;
              }
          }'	

    >>>
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem+"GB"
        maxRetries: 1
    }
    output {
      Array[File] batched_split = glob("*_sample_list_*.txt")
    }
}

task run_split {
    input {
        File combined_vcf
        Int mem
        File batch_file 
    }

    command <<<
      set -e
      ulimit -n 800000
      #para_split() {
      #    input_vcf=$1 
      #    batch_file=$2
      #    out_dir=${batch_file%.txt}
      #    mkdir -p ${out_dir}
      #    echo "Start split ${batch_file} into ${out_dir}"
      #    bcftools plugin split ${input_vcf} -G ${batch_file} -Oz -o ${out_dir}/ 
      #    echo "Start index splitted files for ${batch_file}..."
      #    for f in ${out_dir}/*.vcf.gz; do tabix -p vcf $f; done
      #    echo "${batch_file} is done"
      #}
      

     echo "start split"    
     mkdir ./split_by_samples_hg38
     bcftools plugin split ~{combined_vcf} -G ~{batch_file} -Oz -o ./split_by_samples_hg38/
     echo "start index" 
     for f in ./split_by_samples_hg38/*.vcf.gz; do tabix -p vcf $f; done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem+"GB"
        disk:"local-disk 100 SSD"
        maxRetries: 1
    }

    output {
      Array[File] batched_vcfs = glob("./*/*.vcf.gz")
      Array[File] batched_vcf_idxs = glob("./*/*.vcf.gz.tbi")
    }

}
