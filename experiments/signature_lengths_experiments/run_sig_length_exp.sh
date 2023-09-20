#!/bin/bash
# exit when any command fails
set -e

DATASET="$1"
SIG_LENGTH="$2"

# Define the arrays
scoring_methods=('adjusted_neighborhood_scoring' 'seurat_scoring' 'seurat_ag_scoring' 'seurat_lvg_scoring' 'scanpy_scoring' 'jasmine_scoring_lh' 'jasmine_scoring_or' 'ucell_scoring')
# Loop over scoring methods
for method in "${scoring_methods[@]}"
do
    echo "Starting experiment with dataset=$DATASET, norm_method=mean, and scoring_method=$method ..."
    python signature_length_experiment.py --dataset "$DATASET" --norm_method mean --dge_on_all pseudobulk --scoring_method "$method" --max_sig_length "$SIG_LENGTH"
done

