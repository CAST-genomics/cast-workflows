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
      Array[File] sample_list = run_split.sample_list
      
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
      para_split() {
          input_vcf=$1 
          batch_file=$2
          out_dir=${batch_file%.txt}
          mkdir -p ${out_dir}
          echo "Start split ${batch_file} into ${out_dir}"
          bcftools plugin split ${input_vcf} -G ${batch_file} -Oz -o ${out_dir}/ 
          echo "Start index splitted files for ${batch_file}..."
          for f in ${out_dir}/*.vcf.gz; do tabix -p vcf $f; done
          echo "${batch_file} is done"
      }
      
			echo "split sample files to batches"    
			n_core=$(nproc)
			bcftools query -l ~{combined_vcf} | \
				awk -v group_size=~{sample_size} -v prefix=~{outprefix} -v cpu_n=${n_core} '
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
						batches_per_file = int((total_batches + (cpu_n-1)) / cpu_n);  # Distribute batches across 8 files

						for (i = 0; i < line_count; i++) {
								batch = int(i / group_size) + 1;  # Compute batch number
								file_idx = int((batch - 1) / batches_per_file) + 1;  # Assign batch to one of 8 files
								file = prefix"_sample_list_" file_idx ".txt";
								print lines[i], "-", prefix"_batch_" batch > file;
						}
				}'   	
      export -f para_split
			bash -c 'for sample_batch in $(ls -d ./*_sample_list_*.txt); do para_split ~{combined_vcf} ${sample_batch} & done; wait'
      echo "All done."

     # echo "extract samples first"
     # # split by batches
     # bcftools query -l ~{combined_vcf} | \
     #   awk -v group_size=~{sample_size} -v prefix=~{outprefix} 'BEGIN {FS=OFS="\t"} {print $1,"-",prefix"_batch_"int(NR / group_size)+1}' > sample_list.txt
     # echo "start split"    
     # mkdir ./split_by_samples_hg38
     # bcftools plugin split ~{combined_vcf} -G sample_list.txt -Oz -o ./split_by_samples_hg38/
     # echo "start index" 
     # for f in ./split_by_samples_hg38/*.vcf.gz; do tabix -p vcf $f; done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-plink2:latest"
        memory: mem+"GB"
        disk:"local-disk 200 SSD"
        maxRetries: 1
    }

    output {
      Array[File] sample_list = glob("*sample_list*.txt")
      Array[File] batched_vcfs = glob("./*/*.vcf.gz")
      Array[File] batched_vcf_idxs = glob("./*/*.vcf.gz.tbi")
    }

}
