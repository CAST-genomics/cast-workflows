# Download the required files from Dropbox
wget -O hg38_corrected.psam "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"
wget -O chr1_hg38.pgen.zst "https://www.dropbox.com/s/k5tul2pe5wytjlz/chr1_hg38.pgen.zst?dl=1"
wget -O chr1_hg38.pvar.zst "https://www.dropbox.com/scl/fi/he43aix3klg1uc5745ggm/chr1_hg38_rs.pvar.zst?rlkey=pkk1a30bhz7qok9meilan1kii&dl=1"


plink2 --zst-decompress chr1_hg38.pgen.zst chr1_hg38.pgen
plink2 --zst-decompress chr1_hg38.pvar.zst chr1_hg38.pvar

# Copy sample file
cp hg38_corrected.psam chr1_hg38.psam

# Convert PLINK2 files to VCF format while keeping only selected samples
plink2 --pfile chr1_hg38 --keep AFR-EAS-EUR-SAS-1kgp-ref-list.txt --export vcf --out AFR-EAS-EUR-SAS-1kgp_chr1

# Compress VCF and create index
bgzip AFR-EAS-EUR-SAS-1kgp_chr1.vcf
tabix -p vcf AFR-EAS-EUR-SAS-1kgp_chr1.vcf.gz
echo "VCF processing complete: AFR-EAS-EUR-SAS-1kgp_chr1.vcf.gz"


# Download the required files from Dropbox
wget -O hgdp.psam "https://www.dropbox.com/s/0zg57558fqpj3w1/hgdp.psam?dl=1"
wget -O hgdp_chr1.pgen.zst "https://www.dropbox.com/s/rn7porgrs0yyc7i/hgdp_chr1.pgen.zst?dl=1"
wget -O hgdp_chr1.pvar.zst "https://www.dropbox.com/s/08zgboud52s9ns1/hgdp_chr1.pvar.zst?dl=1"

plink2 --zst-decompress hgdp_chr1.pgen.zst hgdp_chr1.pgen
plink2 --zst-decompress hgdp_chr1.pvar.zst hgdp_chr1.pvar

cp hgdp.psam hgdp_chr1.psam

plink2 --pfile hgdp_chr1 --keep hgdp-AMERICA-MIDDLE_EAST-ref-list.txt --export vcf --out hgdp-AMERICA-MIDDLE_EAST_chr1

bgzip hgdp-AMERICA-MIDDLE_EAST_chr1.vcf
tabix -p vcf hgdp-AMERICA-MIDDLE_EAST_chr1.vcf.gz
echo "VCF processing complete: hgdp-AMERICA-MIDDLE_EAST_chr1.vcf.gz"
##################################################################################

bcftools merge -m all -Oz -o hgdp-1kgp-ref_chr1.vcf.gz hgdp-AMERICA-MIDDLE_EAST_chr1.vcf.gz AFR-EAS-EUR-SAS-1kgp_chr1.vcf.gz
bcftools index --tbi hgdp-1kgp-ref_chr1.vcf.gz

# wget https://faculty.washington.edu/browning/beagle/beagle.17Dec24.224.jar
# wget https://raw.githubusercontent.com/browning-lab/ukb-phasing/refs/heads/master/maps/plink.chr1.GRCh38.map

##################################################################################

bcftools norm -m -any -Oz -o hgdp-1kgp-ref_chr1_norm.vcf.gz hgdp-1kgp-ref_chr1.vcf.gz
bcftools index --tbi hgdp-1kgp-ref_chr1_norm.vcf.gz

bcftools view -m2 -M2 -v snps hgdp-1kgp-ref_chr1_norm.vcf.gz -Oz -o hgdp-1kgp-ref_chr1_norm.biallelic.vcf.gz
bcftools index --tbi hgdp-1kgp-ref_chr1_norm.biallelic.vcf.gz
bcftools norm -d all -Oz -o hgdp-1kgp-ref_chr1_dedup.biallelic.vcf.gz hgdp-1kgp-ref_chr1_norm.biallelic.vcf.gz
bcftools index --tbi hgdp-1kgp-ref_chr1_dedup.biallelic.vcf.gz

# bcftools query -f '%CHROM\t%POS\t%ID\n' hgdp-1kgp-ref_chr1_dedup.biallelic.vcf.gz | sort | uniq -d
# if return nothing, then proceed to next step

##################################################################################

# export LD_LIBRARY_PATH=/home/jupyter/miniconda/pkgs/libarchive-3.7.4-hfab0078_0/lib:$LD_LIBRARY_PATH
# export CONDA_PKGS_DIRS=$HOME/.conda_pkgs/
# mkdir -p $CONDA_PKGS_DIRS

# Create a new Conda environment with OpenJDK 21
# conda create -n java21_env openjdk=21 -y

# Activate the environment
# conda activate java21_env
##################################################################################
sed -i 's/^chr1/1/' plink.chr1.GRCh38.map

java -Xmx600g -jar beagle.17Dec24.224.jar \
    gt=hgdp-1kgp-ref_chr1_dedup.biallelic.vcf.gz \
    map=plink.chr1.GRCh38.map \
    out=hgdp-1kgp-ref_chr1_biallelic.phased \
    nthreads=96

bcftools index --tbi hgdp-1kgp-ref_chr1_biallelic.phased.vcf.gz

##################################################################################
##################################################################################
#!/bin/bash

# Set variables 
VCF_FILE="acaf_threshold_unrelated.chr1.biallelic.phased.vcf.gz"
OUTPUT_PREFIX="flare-output/flare_output.acaf-unrelated-all.chr1"
CHUNK_SIZE=12000
OUTPUT_DIR="flare-output"

# Create output directory
mkdir -p $OUTPUT_DIR

# Get the total number of samples
TOTAL_SAMPLES=$(bcftools query -l $VCF_FILE | wc -l)

# Function to process each chunk
process_chunk() {
    local start=$1
    local end=$2
    local chunk_id=$3

    # Create a temporary file for the sample IDs of this chunk
    tmp_sample_file=$(mktemp)
    bcftools query -l $VCF_FILE | sed -n "${start},${end}p" > $tmp_sample_file

    # Create separate VCF files for each chunk of sample IDs
    bcftools view -S $tmp_sample_file $VCF_FILE -Oz -o ${OUTPUT_DIR}/chunk_${chunk_id}.vcf.gz

    # Clean up temporary files
    rm $tmp_sample_file
}

# Loop to process each chunk based on the total number of samples and chunk size
chunk_id=0
for (( start=1; start<=$TOTAL_SAMPLES; start+=$CHUNK_SIZE )); do
    end=$(( start + CHUNK_SIZE - 1 ))
    if [ $end -gt $TOTAL_SAMPLES ]; then
        end=$TOTAL_SAMPLES
    fi
    process_chunk $start $end $chunk_id
    ((chunk_id++))
done

echo "Splitting chunks complete."

##################################################################################
##################################################################################
zcat hgdp-1kgp-ref_chr1_biallelic.phased.vcf.gz | sed 's/^1\t/chr1\t/' | bgzip -c > hgdp-1kgp-ref_chr1_biallelic.phased.add_chr.vcf.gz
bcftools index --tbi hgdp-1kgp-ref_chr1_biallelic.phased.add_chr.vcf.gz


sed -i 's/^1 /chr1 /' plink.chr1.GRCh38.map

##################################################################################

# Set variables
JAVA_OPTS="-Xmx512g"
FLARE_JAR="flare.jar"
REF_VCF_FILE="hgdp-1kgp-ref_chr1_biallelic.phased.add_chr.vcf.gz"
REF_PANEL="ucsd-ref.panel"
MAP_FILE="plink.chr1.GRCh38.map"
OUTPUT_PREFIX="flare-output/flare_output.6-ancestry-acaf-unrelated-all.chr1"
OUTPUT_DIR="flare-output"
MERGED_LOG_FILE="${OUTPUT_PREFIX}.merged.log"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to run FLARE on each chunk
process_chunk_with_flare() {
    # Process each chunk in numerical order
    find ${OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | while read chunk_file; do
        chunk_id=$(basename $chunk_file | sed 's/chunk_\([0-9]*\).vcf.gz/\1/')
        
        # Define log file for each chunk
        LOG_FILE="${OUTPUT_PREFIX}.chunk_${chunk_id}.log"
        
        # Run FLARE for the chunk and capture log
        java $JAVA_OPTS -jar $FLARE_JAR \
            ref=$REF_VCF_FILE \
            ref-panel=$REF_PANEL \
            gt=$chunk_file \
            map=$MAP_FILE \
            probs=true \
            em=true \
            out=${OUTPUT_PREFIX}.chunk_${chunk_id} &> "$LOG_FILE"
    done

    echo "FLARE processing complete." 
}

# Function to index and merge VCF files
index_and_merge_vcf_files() {
    # Step 1: Index all VCF files if not already indexed
    echo "Indexing VCF files..."
    for vcf in ${OUTPUT_PREFIX}.chunk_*.anc.vcf.gz; do
        if [ ! -f "$vcf.tbi" ]; then
            echo "Indexing $vcf..."
            bcftools index -t $vcf
        fi
    done

    # Step 2: Merge local ancestry VCF files
    echo "Merging local ancestry VCF files..."
    bcftools merge --file-list <(ls ${OUTPUT_PREFIX}.chunk_*.anc.vcf.gz | sort -V) -O z -o ${OUTPUT_PREFIX}.merged.anc.vcf.gz

    # Step 3: Concatenate and index global ancestry files
    echo "Merging and indexing global ancestry files..."
    zcat $(ls ${OUTPUT_PREFIX}.chunk_*.global.anc.gz | sort -V) | gzip > ${OUTPUT_PREFIX}.merged.global.anc.gz

    echo "Processing and merging complete."
}

# Function to merge all log files
merge_log_files() {
    echo "Merging all log files into ${MERGED_LOG_FILE}..."
    cat ${OUTPUT_PREFIX}.chunk_*.log > "$MERGED_LOG_FILE"
    echo "Log merging complete."
}

# Function to index global ancestry files
index_global_ancestry_files() {
    echo "Indexing global ancestry files..."
    for global_file in ${OUTPUT_PREFIX}.chunk_*.global.anc.gz; do
        if [ ! -f "$global_file.tbi" ]; then
            echo "Indexing $global_file..."
            tabix -s1 -b2 -e2 $global_file
        fi
    done

    # Index merged global ancestry file
    if [ -f "${OUTPUT_PREFIX}.merged.global.anc.gz" ]; then
        echo "Indexing ${OUTPUT_PREFIX}.merged.global.anc.gz..."
        tabix -s1 -b2 -e2 ${OUTPUT_PREFIX}.merged.global.anc.gz
    fi
}

# Function to merge all model files
merge_model_files() {
    echo "Merging all model files into ${MERGED_MODEL_FILE}..."
    cat ${OUTPUT_PREFIX}.chunk_*.model > "$MERGED_MODEL_FILE"
    echo "Model merging complete."
}

# Run the functions in order
process_chunk_with_flare
index_and_merge_vcf_files
merge_log_files
merge_model_files
index_global_ancestry_files

##################################################################################
# Reorder the global ancestry label order
zcat flare_output.6-ancestry-acaf-unrelated-all.chr1.merged.global.anc.gz | awk '
NR==1 {print $1, $7, $3, $5, $4, $2, $6; next} 
{print $1, $7, $3, $5, $4, $2, $6}' OFS="\t" | bgzip > flare_output.6-ancestry-acaf-unrelated-all.chr1.merged.global.anc.reordered.gz

bcftools index --tbi flare_output.6-ancestry-acaf-unrelated-all.chr1.merged.anc.vcf.gz
