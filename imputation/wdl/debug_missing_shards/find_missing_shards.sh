#!/bin/bash

workflow="25599b82-617f-4b0c-bce4-7360bbf9873c"
#gsutil ls $WORKSPACE_BUCKET/cromwell-execution/imputation_batch/${workflow}/call-imputation/shard-*/imputation/*/call-sort_index_beagle/cacheCopy/batch_set1_output.sorted_TR.vcf.gz > list_items.txt
gsutil ls $WORKSPACE_BUCKET/cromwell-execution/imputation_batch/${workflow}/call-imputation/shard-*/imputation/*/call-sort_index_beagle/batch_set1_output.sorted_TR.vcf.gz > list_items.txt

cat list_items.txt | awk -F/ '{print ($8)}' | awk -F- '{print($2)}' | sort -n > list_items_shard_number.txt

seq 0 225 > full_list.txt

# Then we can find the shards with missing items by running
#vimdiff full_list.txt list_items_shard_number.txt

