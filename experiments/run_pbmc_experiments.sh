#!/bin/bash
# exit when any command fails
set -e

PATH_TO_REPO="[FILL HERE PATH TO ANS_SUPPLEMENTARY_INFORMATION REPO]"
BASE_STORING_PATH="[FILL HERE PATH TO WHERE EXPERIMENT LOGS/RESULTS SHOULD BE STORED]"

echo "STARTING CONTROL BIAS AND COMPARABLE SCORE RANGE EXPERIMENTS"
echo ""
#####################################################################################################
# CONTROL BIAS EXPERIMENT
#####################################################################################################
cd "$PATH_TO_REPO/experiments/control_bias"
echo "> Starting control bias experiment"

bash run_control_bias.sh --base_storing_path "$BASE_STORING_PATH/control_bias"|& tee -a "$BASE_STORING_PATH/control_bias/logs_control_bias.txt"

echo "> Finished control bias experiment"
echo ""

######################################################################################################
## COMPARABLE SCORE RANGE EXPERIMENT OLD
######################################################################################################
#cd "$PATH_TO_REPO/experiments/comparable_score_ranges"
#echo "> Starting comparable score range experiments EASY TASK"
#echo "> !!!! NOTE: CHANGE STORING PATH IN experiments/run_exp_b_mono_nk.sh !!!!"
#
#bash run_exp_b_mono_nk.sh |& tee -a "$BASE_STORING_PATH/comparable_score_ranges/logs/logs_easy_task.txt"
#
#echo "> Finished comparable score range experiments EASY TASK"
#echo ""
#
#echo "> Starting comparable score range experiments HARD TASK"
#echo "> !!!! NOTE: CHANGE STORING PATH IN experiments/run_exp_b_subtypes.sh !!!!"
#
#bash run_exp_b_subtypes.sh |& tee -a "$BASE_STORING_PATH/comparable_score_ranges/logs/logs_HARD_task.txt"
#
#echo "> Finished comparable score range experiments HARD TASK"
#echo ""


echo "FINISHED CONTROL BIAS AND COMPARABLE SCORE RANGE EXPERIMENTS"