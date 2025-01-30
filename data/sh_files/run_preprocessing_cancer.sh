#!/bin/bash
# exit when any command fails
set -e

cd ..

NORM_METHOD='mean'
# NORM_METHOD='CP10k'

datasets=('crc' 'escc' 'luad_xing' 'breast_small' 'breast_malignant' \
          'luad_kim_malignant' 'luad_kim_malignant_2' \
          'skin_malignant' 'skin_malignant_2' \
          'ovarian_malignant' 'ovarian_malignant_2' )

for dataset in 'ovarian_malignant_2'; do
    echo "STARTING preprocessing $dataset with norm_method=$NORM_METHOD and sample_based=False."
    python preprocess_cancer_data.py --dataset "$dataset" --norm_method "$NORM_METHOD"
    echo "FINISHED preprocessing."
#    echo ""
#
#    echo "STARTING preprocessing $dataset with norm_method=$NORM_METHOD and sample_based=True."
#    python preprocess_cancer_data.py --dataset "$dataset" --norm_method "$NORM_METHOD" --sample_based
#    echo "FINISHED preprocessing."
#    echo ""
done