#!/bin/bash
# exit when any command fails
set -e

echo "Starting experiment with dataset=crc and scoring_method=jasmine_scoring_or ..."
python signature_length_experiment.py --dataset crc --norm_method mean --dge_on_all pseudobulk --scoring_method jasmine_scoring_or --max_sig_length 550
echo " "

echo "Starting experiment with dataset=escc and scoring_method=jasmine_scoring_or ..."
python signature_length_experiment.py --dataset escc --norm_method mean --dge_on_all pseudobulk --scoring_method jasmine_scoring_or --max_sig_length 450
echo " "

echo "Starting experiment with dataset=luad and scoring_method=jasmine_scoring_or ..."
python signature_length_experiment.py --dataset luad --norm_method mean --dge_on_all pseudobulk --scoring_method jasmine_scoring_or --max_sig_length 388
echo "FINISHED!"