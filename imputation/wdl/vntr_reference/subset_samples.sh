#!/bin/bash

# for 70% train and 30% test:
#num_train=650
#num_test=275

# for 80% train and 20% test:
num_train=740
num_train=185
cat vntr_samples_clean_sorted.txt | shuf | head -n $num_train > vntr_samples_${num_train}_train.txt
cat vntr_samples_clean_sorted.txt | grep -v -f vntr_samples_${num_train}_train.txt > vntr_samples_${num_test}_test.txt

# test the samples are split correctly
cat vntr_samples_${num_train}_train.txt  vntr_samples_${num_test}_test.txt | sort --uniq > test.txt

echo "Num lines in the merged file $(wc -l test.txt) (should be 925)"
echo "difference in original and reconstructed sample size (should be 0) : $(comm -3 test.txt vntr_samples_clean_sorted.txt | wc -l)"
