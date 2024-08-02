#!/bin/bash

workflow="a5d62008-af2c-4474-b8c0-d86f3b068d73"
workflow="6a859277-68af-4f8b-a497-a5f5833b4fb5"
gsutil ls gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/cromwell-execution/imputation_batch/${workflow}/call-imputation/shard-*/imputation/*/call-sort_index_beagle/cacheCopy/batch_set1_output.sorted_TR.vcf.gz > list_items.txt
gsutil ls gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/cromwell-execution/imputation_batch/${workflow}/call-imputation/shard-*/imputation/*/call-sort_index_beagle/batch_set1_output.sorted_TR.vcf.gz >> list_items.txt

cat list_items.txt | awk -F/ '{print ($8)}' | awk -F- '{print($2)}' | sort -n > list_items_shard_number.txt

seq 0 225 > full_list.txt

# Then we can find the shards with missing items by running
#vimdiff full_list.txt list_items_shard_number.txt

