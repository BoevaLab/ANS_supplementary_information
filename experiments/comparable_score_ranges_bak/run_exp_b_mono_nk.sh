#!/bin/bash
# exit when any command fails
set -e

dataset="pbmc_b_mono_nk"
### TODO: CHANGE BASE STORING PATH
base_storing_path="...../experiments/comparable_score_ranges/B_NK_Mono/performance_tbls"
scale_stride=0.05

# Compute comparable score range experiment for B, Mono, and NK cells 
python comparable_score_range_experiment.py --dataset "$dataset" --base_storing_path "$base_storing_path" --nr_sigs 3 --scale_stride "$scale_stride"
python comparable_score_range_experiment.py --dataset "$dataset" --base_storing_path "$base_storing_path" --nr_sigs 3 --overlapping_sigs --scale_stride "$scale_stride"
python comparable_score_range_experiment.py --dataset "$dataset" --base_storing_path "$base_storing_path" --nr_sigs 2 --save_adata
python comparable_score_range_experiment.py --dataset "$dataset" --base_storing_path "$base_storing_path" --nr_sigs 2 --overlapping_sigs --save_adata
