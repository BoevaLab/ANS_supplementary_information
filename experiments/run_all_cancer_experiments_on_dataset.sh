#!/bin/bash
# exit when any command fails
set -e

#DATASET=('crc' 'escc' 'luad_xing' 'luad_atlas' 'melanoma')
DATASET="$1"
SCORING_METHODS=('adjusted_neighborhood_scoring' 'seurat_scoring' 'seurat_ag_scoring' 'seurat_lvg_scoring' 'scanpy_scoring' 'jasmine_scoring_lh' 'jasmine_scoring_or' 'ucell_scoring')
MAX_SIG_LENGTH=20
SIG_LENGTH_NOISE=100
PATH_TO_REPO="[FILL HERE PATH TO ANS_SUPPLEMENTARY_INFORMATION REPO]"
BASE_STORING_PATH="[FILL HERE PATH TO WHERE EXPERIMENT LOGS/RESULTS SHOULD BE STORED]"

echo "STARTING ALL EXPERIMENTS ON $DATASET"
echo ""
#####################################################################################################
# Data composition experiment
#####################################################################################################
cd "$PATH_TO_REPO/experiments/data_composition_experiments"
echo "> Starting Data composition experiment on $DATASET"

bash run_data_comp_exp.sh "$DATASET" |& tee -a "$BASE_STORING_PATH/data_composition_experiments/logs_$DATASET.txt"

echo "> Finished Data composition experiment on $DATASET"
echo ""

#####################################################################################################
# Signature length experiment
#####################################################################################################
cd "$PATH_TO_REPO/experiments/signature_lengths_experiments"
echo "> Starting Signature length experiment on $DATASET with maximum signture length $MAX_SIG_LENGTH"

bash run_sig_length_exp.sh "$DATASET" "$MAX_SIG_LENGTH" |& tee -a "$BASE_STORING_PATH/signature_lengths_experiments/logs_$DATASET.txt"

echo "> Finished Signature length experiment on $DATASET"
echo ""

#####################################################################################################
# Signature noise experiment
#####################################################################################################
cd "$PATH_TO_REPOn/experiments/signature_noise_addition_experiments"
echo "> Starting signature noise experiment on $DATASET with signture length $SIG_LENGTH_NOISE"

bash run_sig_noise_exp.sh "$DATASET" "$SIG_LENGTH_NOISE" |& tee -a "$BASE_STORING_PATH/signature_noise_addition_experiments/logs_$DATASET.txt"

echo "> Finished signature noise experiment on $DATASET"
echo ""
#####################################################################################################
echo "FINISHED ALL EXPERIMENTS ON $DATASET"