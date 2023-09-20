#!/bin/bash
# exit when any command fails
set -e

DATASET="$1"

echo "Running data composition experiment on $DATASET"

python data_composition.py --dataset "$DATASET" --norm_method mean --dge_on_all pseudobulk --save

echo "FINISHED!"
