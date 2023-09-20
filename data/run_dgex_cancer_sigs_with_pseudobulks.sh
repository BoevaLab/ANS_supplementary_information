#!/bin/bash
# exit when any command fails
set -e

NORM_METHOD='mean'
# NORM_METHOD='CP10k'

# for dataset in 'crc' 'escc' 'luad_xing' 'luad_atlas' 'melanoma'; do
# for dataset in 'luad_xing' 'luad_atlas' 'melanoma'; do
 for dataset in 'luad_xing'; do
    echo "STARTING extraction malignant signature for $dataset with norm_method=$NORM_METHOD."
    python DGEX_cancer_sig_with_pseudobulk_and_deseq2.py --dataset $dataset --norm_method "$NORM_METHOD" --save_fig
    echo "FINISHED experiment."
    echo ""
done
