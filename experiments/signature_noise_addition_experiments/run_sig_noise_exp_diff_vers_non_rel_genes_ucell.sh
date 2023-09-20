#!/bin/bash
# exit when any command fails
set -e

max_abs_log2fc=0.5
min_pval=0.01

# Loop over remaining attributes
for dataset in "crc" "escc" "luad"
do
    echo "Running noise experiment with noise genes with |log2FC|<0.5 for dataset=$dataset"
    python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method ucell_scoring  --signature_length 100 --nr_sims 20 --max_abs_log2fc "$max_abs_log2fc"  --min_pval "$min_pval"
done


