#!/bin/bash
# exit when any command fails
set -e

#for dataset in 'crc' 'escc' 'luad_xing' 'luad_atlas' 'melanoma'; do
for dataset in 'luad_xing' 'luad_atlas' 'melanoma'; do
#  for norm_method in 'mean' 'CP10k'; do
    echo "STARTING extraction malignant signature for $dataset with norm_method=mean."
    python DGEX_cancer_sig_on_individual_samples.py --dataset $dataset --norm_method mean
    echo "FINISHED experiment."
#  done
done
