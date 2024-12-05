#!/bin/bash

datasets=("breast_malignant" "luad_kim_malignant" "luad_kim_malignant_2" "skin_malignant" "skin_malignant_2" "ovarian_malignant" "pbmc_b_mono_nk" "pbmc_b_subtypes" "pbmc_cd4_subtypes" "pbmc_cd8_subtypes")
gt_cols=("gene_module" "Cell_subtype" "Cell_subtype" "level2_celltype" "level2_celltype" "cluster_label" "celltype.l1" "celltype.l2" "celltype.l2" "celltype.l2")
sample_cols=("Patient" "sample" "sample" "patient" "patient" "patient_id" "orig.ident" "orig.ident" "orig.ident" "orig.ident")

# Create logs directory if it doesn't exist
mkdir -p "./logs"

for i in "${!datasets[@]}"; do
    dataset="${datasets[$i]}"
    gt_col="${gt_cols[$i]}"
    sample_col="${sample_cols[$i]}"

    echo "Processing dataset: $dataset with gt_col: $gt_col and sample_col: $sample_col"

    # Base method
    echo "Running base method..."
    python comparable_score_range_experiment.py \
       --dataset "${dataset}" \
       --sample_col "${sample_col}" \
       --gt_annotation_col "${gt_col}" \
       --save_signatures \
       --verbose > "logs/${dataset}_base.log" 2>&1

    # Remove overlapping genes
    echo "Running with remove_overlapping_genes..."
    python comparable_score_range_experiment.py \
           --dataset "${dataset}" \
           --sample_col "${sample_col}" \
           --gt_annotation_col "${gt_col}" \
           --remove_overlapping_genes \
           --verbose > "logs/${dataset}_remove_overlapping.log" 2>&1

    # Use gene pool
    echo "Running with use_gene_pool..."
    python comparable_score_range_experiment.py \
           --dataset "${dataset}" \
           --sample_col "${sample_col}" \
           --gt_annotation_col "${gt_col}" \
           --use_gene_pool \
           --verbose > "logs/${dataset}_gene_pool.log" 2>&1

    # Both modifications
    echo "Running with both modifications..."
    python comparable_score_range_experiment.py \
           --dataset "${dataset}" \
           --sample_col "${sample_col}" \
           --gt_annotation_col "${gt_col}" \
           --remove_overlapping_genes \
           --use_gene_pool \
           --verbose > "logs/${dataset}_both_mods.log" 2>&1
done