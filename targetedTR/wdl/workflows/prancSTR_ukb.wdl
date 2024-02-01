version 1.0

#### DNANexus
# WDL script for UKBiobank data processing using prancSTR
# Currently aiming to do analysis on all STR's genome wide for a few samples

import "../tasks/dumpstr.wdl" as dumpstr_t
import "../tasks/prancstr.wdl" as prancstr_w

workflow prancSTR {
	input {
		File input_vcf
        String vcftype
        String readfield
        String out_prefix
		String region
	}

	### Running prancSTR on final hipstr vcf ###
	call prancstr_w.run_prancstr as prancstr{
		input:
			input_vcf=input_vcf,
			vcftype=vcftype,
			readfield=readfield,
			region=region,
			out_prefix=out_prefix
	}

	### Output files ####
	output {
		File tabvcf = prancstr.outfile
	}
}
