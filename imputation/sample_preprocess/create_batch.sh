'''example code
gsutil cp $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample/aou_sample_list.txt .
./create_batch.sh aou_sample_list.txt
gsutil cp -r sample_batch* $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample_batch/
#generate sample_manifest.csv
gsutil ls $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample_batch/ > sample_manifest.txt
sed 's/$/,/' sample_manifest.txt  > sample_manifest.csv
#save in google cloud
gsutil cp sample_manifest.csv $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample/
'''

SAMPLEFILE=$1

split -l 1000 $1 sample_batch
counter=1
for file in sample_batch*
do
    mv "$file" "sample_batch$counter"
    ((counter++))
done