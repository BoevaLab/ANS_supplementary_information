#!/bin/bash
# exit when any command fails
set -e

DATASET="$1"
SIG_LENGTH="$2"

max_abs_log2fc=0.5
min_pval=0.01

# Define the arrays
scoring_methods=('adjusted_neighborhood_scoring' 'seurat_scoring' 'seurat_ag_scoring' 'seurat_lvg_scoring' 'scanpy_scoring' 'jasmine_scoring_lh' 'jasmine_scoring_or' 'ucell_scoring')

# Loop over scoring methods
for method in "${scoring_methods[@]}"
do
    echo "> Starting experiment with dataset=$DATASET, scoring_method=$method, and  signature_length_1=$SIG_LENGTH ..."
    # echo ">> all remaining genes as noise"
    # python signature_noise_experiment.py --dataset "$DATASET" --norm_method mean --dge_on_all pseudobulk --scoring_method "$method" --signature_length "$SIG_LENGTH" --nr_sims 20
    echo ">> genes with max_abs_log2fc=$max_abs_log2fc and min_pval=$min_pval as noise"
    python signature_noise_experiment.py --dataset "$DATASET" --norm_method mean --dge_on_all pseudobulk --scoring_method "$method" --signature_length "$SIG_LENGTH" --nr_sims 20 --max_abs_log2fc "$max_abs_log2fc" --min_pval "$min_pval"
    # echo ">> genes with min_pval=$min_pval as noise"
    # python signature_noise_experiment.py --dataset "$DATASET" --norm_method mean --dge_on_all pseudobulk --scoring_method "$method" --signature_length "$SIG_LENGTH" --nr_sims 20 --min_pval "$min_pval"
    echo ">> FINISHED with scoring_method=$method!"
    echo " "
done

