#!/bin/bash

# Deprecated
#wget "https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"

# Download
#wget "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"

# Reformat
for chr_idx in $(seq 1 22) "X"; do
	cat "plink.chr${chr_idx}.GRCh38.map" | awk '{print("chr" $0)}' > "plink.chr${chr_idx}_w_chr.GRCh38.map"
done
