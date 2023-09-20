#!/bin/bash
# exit when any command fails
set -e

NORM_METHOD='mean'
# NORM_METHOD='CP10k'

#for dataset in 'crc' 'escc' 'luad_xing' 'luad_atlas' 'melanoma'; do
#for dataset in 'luad_xing' 'melanoma'; do
for dataset in 'breast_large' 'breast_small'; do
    echo "STARTING preprocessing $dataset with norm_method=$NORM_METHOD and sample_based=False."
    python preprocess_cancer_data.py --dataset "$dataset" --norm_method "$NORM_METHOD"
    echo "FINISHED preprocessing."

    echo ""

    echo "STARTING preprocessing $dataset with norm_method=$NORM_METHOD and sample_based=True."
    python preprocess_cancer_data.py --dataset "$dataset" --norm_method "$NORM_METHOD" --sample_based
    echo "FINISHED preprocessing."
    echo ""
done