#!/bin/bash

vcf_file=$1
bcftools view -i 'ID="."' $vcf_file | wc -l
