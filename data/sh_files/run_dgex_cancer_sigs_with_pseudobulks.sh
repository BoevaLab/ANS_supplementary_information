#!/bin/bash
# exit when any command fails
set -e

NORM_METHOD='mean'
cd ../malignant_cell_signatures

datasets=('crc' 'escc' 'luad_xing')

for dataset in ${datasets[@]}; do
    echo "STARTING extraction malignant signature for $dataset with norm_method=$NORM_METHOD."
    python DGEX_cancer_sig_with_pseudobulk_and_deseq2.py --dataset $dataset --norm_method "$NORM_METHOD" --save_fig
    echo "FINISHED experiment."
    echo ""
done
