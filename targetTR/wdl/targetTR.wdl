version 1.0

import "hipstr_multi.wdl" as hipstr_multi_t
import "merge_hipstr.wdl" as merge_t
import "dumpstr.wdl" as dumpstr_t

workflow targetTR {
	input {
		String outprefix
		Array[Array[File]] cram_file_batches = []
		Array[Array[File]] cram_index_batches = []
		File genome
		File genome_index
		File tr_bed
		Boolean infer_samps_from_file = false
		Array[File] cram_file_batches_str = []
		String GOOGLE_PROJECT = ""
		String GCS_OAUTH_TOKEN = ""
		Float sleep_const = 0
		Boolean separate_hipstr_runs = false
		Boolean longtr = false
		Int? hipstr_mem
		Int? merge_mem
		String? extra_hipstr_args
	}

	### Call HipSTR on batches of samples ###
	Boolean using_batch_files = (length(cram_file_batches_str)>0)
	Int num_batches = if using_batch_files
		then length(cram_file_batches_str)  
		else length(cram_file_batches)
	scatter(i in range(num_batches)) {
		if (using_batch_files) {
			File crams_file = cram_file_batches_str[i]
		}
		if (!using_batch_files) {
			Array[File] crams = cram_file_batches[i]
			Array[File] cram_indices = cram_index_batches[i]
		}
		Int sleep_seconds = ceil(i*sleep_const)
		call hipstr_multi_t.run_hipstr as run_hipstr {
			input :
				bams=crams,
				bam_indices=cram_indices,
				genome=genome,
				genome_index=genome_index,
				str_ref=tr_bed,
				out_prefix=outprefix+".BATCH"+i,
				infer_samps_from_file = infer_samps_from_file,
				using_batch_files = using_batch_files,
				bams_file = crams_file,
				GOOGLE_PROJECT = GOOGLE_PROJECT,
				GCS_OAUTH_TOKEN = GCS_OAUTH_TOKEN,
				sleep_seconds = sleep_seconds,
				hipstr_mem = hipstr_mem,
				separate_hipstr_runs = separate_hipstr_runs,
				longtr = longtr,
				extra_hipstr_args = extra_hipstr_args
		}
	}

	### Merge HipSTR across all samples ###
	call merge_t.merge_hipstr as merge_hipstr {
		input :
			vcfs=run_hipstr.outfile,
			vcf_indexes=run_hipstr.outfile_index,
			out_prefix=outprefix,
			merge_mem=merge_mem
	}

	### DumpSTR on merged VCF ###
	call dumpstr_t.run_dumpstr as dumpstr {
		input :
			vcf=merge_hipstr.outfile,
			out_prefix=outprefix+".filtered"
	}

	### Zip and index the VCF ###
	call sort_index {
		input :
			vcf=dumpstr.outfile
	}

	### Output files ####
	output {
		File finalvcf = sort_index.outvcf
		File finalvcf_index = sort_index.outvcf_index
	}
}

task sort_index {
	input {
		File vcf
	}

	String basename = basename(vcf, ".vcf")

	command <<<
		vcf-sort ~{vcf} | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
	>>>

	runtime {
		docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
	}

	output {
		File outvcf = "${basename}.sorted.vcf.gz"
		File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
	}
}
