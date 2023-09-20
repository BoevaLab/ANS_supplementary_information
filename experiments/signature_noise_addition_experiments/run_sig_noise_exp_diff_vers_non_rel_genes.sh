#!/bin/bash
# exit when any command fails
set -e

# Function to run experiments for a given scoring method
run_experiment() {
    local scoring_method=$1
    local dataset=$2
    local signature_length_1=$3
    local signature_length_2=$4

    echo "> Starting experiment with dataset=$dataset and scoring_method=$scoring_method,  signature_length_1=$signature_length_1,  signature_length_2=$signature_length_2 ..."
    # echo ">> all remaining genes as noise"
    # python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method "$scoring_method" --signature_length "$signature_length_1" --nr_sims 20
    # python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method "$scoring_method" --signature_length "$signature_length_2" --nr_sims 20
    echo ">> genes with max_abs_log2fc=$max_abs_log2fc and min_pval=$min_pval as noise"
    python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method "$scoring_method" --signature_length "$signature_length_1" --nr_sims 20 --max_abs_log2fc "$max_abs_log2fc" --min_pval "$min_pval"
    python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method "$scoring_method" --signature_length "$signature_length_2" --nr_sims 20 --max_abs_log2fc "$max_abs_log2fc" --min_pval "$min_pval"
    # echo ">> genes with min_pval=$min_pval as noise"
    # python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method "$scoring_method" --signature_length "$signature_length_1" --nr_sims 20 --min_pval "$min_pval"
    # python signature_noise_experiment.py --dataset "$dataset" --norm_method mean --dge_on_all pseudobulk --scoring_method "$scoring_method" --signature_length "$signature_length_2" --nr_sims 20 --min_pval "$min_pval"
    echo ">> FINISHED with scoring_method=$scoring_method!"
    echo " "
}

max_abs_log2fc=0.5
min_pval=0.01

echo "Running experiments for the scoring methods: [$@]"

dataset="$1"
signature_length_1="$2"
signature_length_2="$3"
shift 3

# Loop over remaining attributes
for sc_method in "$@"
do
    run_experiment "$sc_method" "$dataset" "$signature_length_1" "$signature_length_2"
done


