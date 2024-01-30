#!/bin/bash

if [ $# != 2 ]; then
	echo "Usage: bash $0 <sample_id> <region>"
	echo "Region is the coordinates for target region e.g. chr1:1-100"
	exit 1
fi
sample_id=$1
region=$2

# Set flags for samtools to read over gcloud URL
export HTSLIB_CONFIGURE_OPTIONS="--enable-gs"
export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"
export GCS_REQUESTER_PAYS_PROJECT="$GOOGLE_PROJECT"

# Alignment file extracted based on sample_id.
alignment_file="gs://fc-aou-datasets-controlled/pooled/longreads/v7_base/bam/${sample_id}/GRCh38/${sample_id}.bam"

# Locally installed samtools updated version >=1.13.
local_samtools="/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/vntr_workspace/packages/bin/bin/samtools"


# Run samtools extracting target region and create corresponding index file.
$local_samtools view -h -b -o target_input_unsorted.bam ${alignment_file} ${region}
$local_samtools sort -o target_input.bam target_input_unsorted.bam
$local_samtools index target_input.bam

# Copy the resulting bam file and index file to the project bucket.
gsutil -u $GOOGLE_PROJECT cp target_input.bam ${WORKSPACE_BUCKET}/saraj/data/target_input.bam
gsutil -u $GOOGLE_PROJECT cp target_input.bam.bai ${WORKSPACE_BUCKET}/saraj/data/target_input.bam.bai
