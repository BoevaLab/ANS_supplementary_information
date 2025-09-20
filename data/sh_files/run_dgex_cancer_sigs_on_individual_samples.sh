#!/bin/bash
# exit when any command fails
set -e

NORM_METHOD='mean'
cd ../malignant_cell_signatures

datasets=('crc' 'escc')

for dataset in ${datasets[@]}; do
    echo "STARTING extraction malignant signature for $dataset with norm_method=mean."
    python DGEX_cancer_sig_on_individual_samples.py --dataset $dataset --norm_method $NORM_METHOD
    echo "FINISHED experiment."
#  done
done
