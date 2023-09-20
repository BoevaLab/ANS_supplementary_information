#!/bin/bash
# exit when any command fails
set -e

# echo "Starting experiment with dataset=crc and scoring_method=$1 ..."
# python signature_noise_experiment.py --dataset crc --norm_method mean --dge_on_all pseudobulk --scoring_method "$1" --signature_length 100 --nr_sims 20
# echo " "

echo "Starting experiment with dataset=escc and scoring_method=$1 ..."
python signature_noise_experiment.py --dataset escc --norm_method mean --dge_on_all pseudobulk --scoring_method "$1" --signature_length 100 --nr_sims 20
echo " "

echo "Starting experiment with dataset=luad and scoring_method=$1 ..."
python signature_noise_experiment.py --dataset luad --norm_method mean --dge_on_all pseudobulk --scoring_method "$1" --signature_length 100 --nr_sims 20
echo "FINISHED!"
