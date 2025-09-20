#!/bin/bash
# exit when any command fails
set -e

cd ../preprocessing

NORM_METHOD='mean'

datasets=(
    'crc' 
    'escc' 
    'luad_xing' 
    'breast_basal_like_samples' 
    'breast_malignant' 
    'luad_kim'
    'luad_kim_malignant' 
    'luad_kim_malignant_3ca' 
    'skin_malignant' 
    'skin_malignant_manual' 
    'ovarian_malignant' 
    'ovarian_malignant_cellxgene' 
)

for dataset in ${datasets[@]}; do
    echo "STARTING preprocessing $dataset with norm_method=$NORM_METHOD and sample_based=False."
    python preprocess_cancer_data.py --dataset "$dataset" --norm_method "$NORM_METHOD"
    echo "FINISHED preprocessing."
done