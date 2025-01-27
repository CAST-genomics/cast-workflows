#!/bin/bash

VCF=$1
OUTPREFIX=$2
STR=$3

# Extract header (chrom, position, ref, alt, sample genotypes)
bcftools view -h $VCF | grep "^#CHROM" | cut -f 1,2,3,4,5,10- > ${OUTPREFIX}.txt

# Extract imputation genotype (variant ID and sample genotypes)
bcftools query -R $STR -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%TGT\t]\n" $VCF >> ${OUTPREFIX}.txt

# Prepare header for output file with sample IDs
echo -e "variant\t ref \t$(bcftools query -R $STR -f"%ID\t" $VCF)" > ${OUTPREFIX}_gtsum.txt

# Extract genotype information (remove header lines, we only want the g
cat ${OUTPREFIX}.txt | cut -f 3,4,6- > ${OUTPREFIX}_cut.txt

# Compute copy number sum per variant for each sample
echo "Calculating copy number sum per variant for each sample..."

awk 'BEGIN {OFS="\t"}
{
    # Print the header line (only for the first record)
    if (NR == 1) {
        print;  # Print the header
    }

    # For the remaining lines (after the header), calculate the copy number sum for each variant and sample
    if (NR > 1) {
        variant = $1;  # The first column is the ID
        ref = length($2);   # The second column is the ref allele length 

        printf "%s\t%s", variant, ref;

        # Loop through each sample genotype (starting from the third column)
        for(i=3; i<=NF; i++) {  
            sum = 0;

            # Handle the case where genotype is in the format of xx|xx
            if ($i ~ /\|/) {
                split($i, alleles, "|");  # Split the genotype by "|"
                sum = length(alleles[1]) + length(alleles[2]);  # Sum the lengths of the two alleles
            } else {
                sum = $i;  
            }

            # Output the copy number sum for this sample
            printf("\t%d", sum);
        }
        print "";  # Print a newline after processing each variant
    }
}' ${OUTPREFIX}_cut.txt > ${OUTPREFIX}_gtsum.txt