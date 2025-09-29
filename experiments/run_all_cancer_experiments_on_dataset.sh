#!/bin/bash
# exit when any command fails
set -e

DATASET="$1"
SCORING_METHODS=('adjusted_neighborhood_scoring' 'seurat_scoring' 'seurat_ag_scoring' 'seurat_lvg_scoring' 'scanpy_scoring' 'jasmine_scoring_lh' 'jasmine_scoring_or' 'ucell_scoring')
MAX_SIG_LENGTH=20
SIG_LENGTH_NOISE=100
# PATH_TO_REPO="[FILL HERE PATH TO ANS_SUPPLEMENTARY_INFORMATION REPO]"
# BASE_STORING_PATH="[FILL HERE PATH TO WHERE EXPERIMENT LOGS/RESULTS SHOULD BE STORED]"
PATH_TO_REPO="/Users/lciernik/Documents/TUB/projects/ans_scoring/ANS_supplementary_information"
BASE_STORING_PATH="/Users/lciernik/Documents/TUB/projects/ans_scoring/reproduce_project/results"

echo "STARTING ALL EXPERIMENTS ON $DATASET"
echo ""
#####################################################################################################
# Data composition experiment
#####################################################################################################
cd "$PATH_TO_REPO/experiments/data_composition_experiments"
echo "> Starting Data composition experiment on $DATASET"
log_path_file="$BASE_STORING_PATH/data_composition_experiments"
if [ ! -d "$log_path_file" ]; then
  mkdir -p "$log_path_file"
fi

bash run_data_comp_exp.sh "$DATASET" 2>&1 | tee -a "$log_path_file/logs_$DATASET.txt"

echo "> Finished Data composition experiment on $DATASET"
echo ""

# #####################################################################################################
# # Signature length experiment
# #####################################################################################################
# cd "$PATH_TO_REPO/experiments/signature_lengths_experiments"
# echo "> Starting Signature length experiment on $DATASET with maximum signture length $MAX_SIG_LENGTH"
# log_path_file="$BASE_STORING_PATH/signature_lengths_experiments"
# if [ ! -d "$log_path_file" ]; then
#   mkdir -p "$log_path_file"
# fi

# bash run_sig_length_exp.sh "$DATASET" "$MAX_SIG_LENGTH" 2>&1 | tee -a "$log_path_file/logs_$DATASET.txt"



# echo "> Finished Signature length experiment on $DATASET"
# echo ""

# #####################################################################################################
# # Signature noise experiment
# #####################################################################################################
# cd "$PATH_TO_REPO/experiments/signature_noise_addition_experiments"
# echo "> Starting signature noise experiment on $DATASET with signture length $SIG_LENGTH_NOISE"
# log_path_file="$BASE_STORING_PATH/signature_noise_addition_experiments"
# if [ ! -d "$log_path_file" ]; then
#   mkdir -p "$log_path_file"
# fi

# bash run_sig_noise_exp.sh "$DATASET" "$SIG_LENGTH_NOISE" 2>&1 | tee -a "$log_path_file/logs_$DATASET.txt"

# echo "> Finished signature noise experiment on $DATASET"
# echo ""
# #####################################################################################################
# echo "FINISHED ALL EXPERIMENTS ON $DATASET"