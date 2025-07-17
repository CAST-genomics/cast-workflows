#!/bin/bash

for chr1 in {22..1}
do
  for chr2 in {21..1}
  do
    if [ "$chr1" -eq "$chr2" ]; then
      continue
    fi
    for c in eur eas nat afr sas mid
    do
      echo "$chr1" "$chr2" "$c" "afr"
      #python make_corr_plot.py "$chr1" "$chr2" "$c" "afr"
    done
  done
done

