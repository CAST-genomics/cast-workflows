# TargetTR on All of Us workbench

All of this must be run on the workbench itself. Once in the workbench, launch a terminal session and to the following.

1. Clone the github and navigate to the launch directory.

```
cd workspaces/impactofglobalandlocalancestryongenomewideassociationv6studies/
mkdir mgymrek-workspace
cd mgymrek-workspace/
git clone https://github.com/cast-genomics/cast-workflows/
cd cast-workflows/
git checkout targetedTR
cd targetedTR/launch
```

2. Set up
```
curl https://github.com/broadinstitute/cromwell/releases/download/77/womtool-77.jar -o womtool-77.jar -L
curl https://github.com/broadinstitute/cromwell/releases/download/77/cromwell-77.jar -o cromwell-77.jar -L
curl -s "https://get.sdkman.io" -o install_sdkman.sh
bash install_sdkman.sh
source "/home/jupyter/.sdkman/bin/sdkman-init.sh"
sdk install java 11.0.14-tem
```

3. Copy the `manifest.csv` file with the CRAM paths to the launch directory.
```
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/cram/manifest.csv .
```

4. Run a test job:

```
./targetTR_launcher_aou.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 2 \
  --batch-num 2 \
  --file-list manifest.csv
```
