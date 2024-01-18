#!/bin/bash
# This script has to be run every time a new environment and disk is allocated.
# This is an adaptive version of targetedTR/launch_aou/README_aou.md based on advntr.

# Clone the cast-workflows repository
cd workspaces/impactofglobalandlocalancestryongenomewideassociationv6studies/
mkdir vntr-workspace; cd vntr-workspace/
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/
git checkout workflow-advntr
cd advntr/wdl/tasks

# Install cromwell v86. This version is needed for "localization_optional" feature referenced in advntr.wdl.
curl https://github.com/broadinstitute/cromwell/releases/download/86/womtool-86.jar -o womtool-86.jar -L
curl https://github.com/broadinstitute/cromwell/releases/download/86/cromwell-86.jar -o cromwell-86.jar -L

# Install other dependencies.
curl -s "https://get.sdkman.io" -o install_sdkman.sh
bash install_sdkman.sh
source "/home/jupyter/.sdkman/bin/sdkman-init.sh"
sdk install java 11.0.14-tem

# Copy the manifest file for long read WGS data.
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/long_read/manifest.csv .
