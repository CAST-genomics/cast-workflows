#!/bin/bash

gsutil ls gs://fc-secure-f6524c24-64d9-446e-8643-415440f52b46/cromwell-execution/imputation_batch/7820fb80-a03f-4544-a9b3-08f2b830a183/call-imputation/shard-*/imputation/*/call-sort_index_beagle/batch_set1_output.sorted.vcf.gz > list_items.txt

cat list_items.txt | awk -F/ '{print ($8)}' | awk -F- '{print($2)}' | sort -n > list_items_shard_number.txt

seq 0 225 > full_list.txt

# Then we can find the shards with missing items by running
#vimdiff full_list.txt list_items_shard_number.txt

