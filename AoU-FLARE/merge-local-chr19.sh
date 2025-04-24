#!/bin/bash

# Special case handling: merge FLARE local ancestry calling files on Chr19
# Standard workflow failed (only on *.anc.vcf.gz files merging on Chr19)

set -euo pipefail

# Output merged VCF
OUTPUT_VCF="flare_output.seperated-ancestry-run-acaf-unrelated.chr19.merged.anc.vcf.gz"

# Temporary directory for fixed files
TEMP_DIR="tmp_merge_final"
mkdir -p "$TEMP_DIR"

# Create a file list of all relevant anc.vcf.gz files
VCF_LIST="$TEMP_DIR/vcf_list.txt"
find ancestry-flare-EUR-run ancestry-flare-AMR-run ancestry-flare-AFR-run -name "*.anc.vcf.gz" | sort -V > "$VCF_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-EAS.chr19.anc.vcf.gz" >> "$VCF_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-MID.chr19.anc.vcf.gz" >> "$VCF_LIST"
echo "flare_output.six-ancestry-acaf-unrelated-SAS.chr19.anc.vcf.gz" >> "$VCF_LIST"
3.

echo "Merging the following VCF files:"
cat "$VCF_LIST"

# Array to hold the paths of "fixed" VCFs
MODIFIED_VCFS=()

# For each original VCF, inject the 'END' header if missing, then compress + index
while read -r vcf; do
    base=$(basename "$vcf" .vcf.gz)
    mod_vcf_gz="$TEMP_DIR/${base}.fixed.vcf.gz"

    (
        # 1) Print header, inserting END definition if absent
        bcftools view -h "$vcf" | \
        awk 'BEGIN { inserted=0 }
             /^##INFO=<ID=END,/ { inserted=1 }
             /^#CHROM/ && inserted==0 {
                 print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">"
                 inserted=1
             }
             { print }'

        # 2) Print the variant lines (no header)
        bcftools view -H "$vcf"
    ) | bgzip -c > "$mod_vcf_gz"

    # Index the newly created fixed VCF
    tabix -f -p vcf "$mod_vcf_gz"

    MODIFIED_VCFS+=("$mod_vcf_gz")
done < "$VCF_LIST"

# Merge the "fixed" VCFs
printf "%s\n" "${MODIFIED_VCFS[@]}" > "$TEMP_DIR/modified_vcf_list.txt"
bcftools merge --file-list "$TEMP_DIR/modified_vcf_list.txt" -Oz -o "$OUTPUT_VCF"

# Index the merged file
bcftools index --tbi "$OUTPUT_VCF"

echo "Final merged VCF created and indexed: $OUTPUT_VCF"
