#!/bin/bash

subtypes=("B memory kappa" "B naive kappa" "B naive lambda" "CD14 Mono" "CD16 Mono" "CD4 CTL" "CD4 Naive" "CD4 TCM_1" "CD4 TCM_3" "CD4 TEM_1" "CD4 TEM_3" "CD8 Naive" "CD8 TEM_1" "CD8 TEM_2" "CD8 TEM_4" "CD8 TEM_5" "MAIT" "NK_1" "NK_2" "NK_3" "Platelet" "cDC2_2")

# Sublist 1
# subtypes=("B memory kappa" "B naive kappa" "B naive lambda" "CD14 Mono" "CD16 Mono" "CD4 CTL" "CD4 Naive")

# Sublist 2
# subtypes=("CD4 TCM_1" "CD4 TCM_3" "CD4 TEM_1" "CD4 TEM_3" "CD8 Naive" "CD8 TEM_1" "CD8 TEM_2")

# Sublist 3
# subtypes=("CD8 TEM_4" "CD8 TEM_5" "MAIT" "NK_1" "NK_2" "NK_3" "Platelet" "cDC2_2")

#conda activate scrna-env

for subtype in "${subtypes[@]}"; do
    echo "Computing for subtype: $subtype"
    python compute_control_bias.py --subtype "$subtype"
    echo ""
    
done
