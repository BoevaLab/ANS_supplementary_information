# Robust and efficient annotation of cell states through gene signature scoring

This repository accompanies the research paper:

**Laure Ciernik\*, Agnieszka Kraft\*, Florian Barkmann, Joséphine Yates, and Valentina Boeva**  
*"Robust and efficient annotation of cell states through gene signature scoring"*  
📄 DOI: [https://doi.org/10.1101/2023.09.20.558114](https://doi.org/10.1101/2023.09.20.558114)

*\*Equal contribution*

## 📋 Table of Contents
- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Data Setup](#data-setup)
- [Repository Structure](#repository-structure)
- [Running Experiments](#running-experiments)
- [Correspondance](#correspondance)

## 📝 Overview

This repository contains Jupyter notebooks and scripts for reproducibility of experiments and visualizations from our research on gene signature scoring methods. 

> **⚠️ Important**: This code is shared for reproducibility purposes and is **not organized as a package**. For the production-ready implementation of ANS, please use our [packaged repository](https://github.com/BoevaLab/ANS_signature_scoring).

## 🚀 Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/BoevaLab/ANS_supplementary_information.git
   cd ANS_supplementary_information
   ```

2. **Set up environment**
   ```bash
   conda create -n scoring_env python=3.9.18
   conda activate scoring_env
   pip install -r requirements.txt
   ```
3. **Set up project folder**
   - Edit `data/constants.py` and replace all `TODO` placeholders with your actual paths
   - Update path placeholders in notebooks and scripts marked with `TODO`
   ```bash
   python data/constants.py
   ```
4. **Download data**
   ```bash
   cd data
   python download_preprocessed_datasets.py
   ```
5. **Install ANS package** (if not already installed)
   ```bash
   pip install git+https://github.com/BoevaLab/ANS_signature_scoring.git
   ```

## 🛠 Installation

### Prerequisites
- Python 3.9.18
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- R (for specific experiments - see requirements below)

### Environment Setup

Choose one of the following methods:

**Option 1: Using pip (recommended)**
```bash
conda create -n scoring_env python=3.9.18
conda activate scoring_env
pip install -r requirements.txt
```

**Option 2: Using conda environment file**
```bash
# With build hashes (exact reproduction)
conda env create -f environment_with_build_hash.yml

# Without build hashes (more flexible)
conda env create -f environment.yml
```

### R Dependencies
For R-based experiments, install the packages listed in `notebooks/construction_scoring_methods/session_info.txt`.

## 📁 Data Setup

### Path Configuration
1. **Set base path**: In `data/constants.py`, set `BASE_PATH_DRIVE` to your project folder
2. **Required structure**:
   ```
   your_project_folder/
   ├── data/           # Downloaded datasets go here
   ├── experiments/    # Experiment outputs
   └── results/        # Results outputs (logs, potentially other results)
   ```
3. **Replace placeholders**: Search for `TODO` in all files and replace with actual paths

### Available Datasets
- CRC, ESCC, LUAD (Xing), LUAD (Kim)
- Breast cancer, skin cancer, ovarian cancer
- PBMC datasets
- CanSig preprocessed datasets

#### Automatic Datasets Download
```bash
python data/download_preprocessed_datasets.py
```
#### Manual Datasets Download
Download all preprocessed datasets and annotations from our [Google Drive folder](https://drive.google.com/drive/folders/10L2gqapJbyOn_MbrZRHQG--n0Xj7wIyg?usp=sharing).


## 🗂 Repository Structure

```
├── README.md
├── data/                              # Data loading and preprocessing
│   ├── constants.py                   # ⚠️ Configure paths here
│   ├── load_data.py
│   ├── download_preprocessed_datasets.py
│   └── run_preprocessing_cancer.sh
├── experiments/                       # Experiment scripts
│   ├── comparable_score_ranges/
│   ├── control_bias/
│   ├── data_composition_experiments/
│   ├── runtime/
│   ├── signature_lengths_experiments/
│   ├── signature_noise_addition_experiments/
│   ├── run_all_experiments_on_dataset.sh
│   └── run_pbmc_experiments.sh
├── notebooks/                         # Analysis notebooks
│   ├── EMT_scoring_experiments/       # 📊 Case study: EMT signal decoupling
│   ├── comparable_score_ranges/
│   ├── construction_scoring_methods/
│   ├── control_genes_selection_experiments/
│   ├── correlation_scores_with_TRC_and_MTP_experiments/
│   ├── data_composition_experiments/
│   ├── signature_lengths_experiments/
│   └── signature_noise_addition_experiments/
├── environment.yml
├── environment_with_build_hash.yml
└── requirements.txt
```
The `experiments` and `notebooks` folders contain the code for executed experiments. While the first contains the python scripts and R scripts, to run the experiments, the notebooks are used to create the visualizations. See "Rerunning experiments" below for details on how to rerun the experiments.

## 🔬 Running Experiments
All experiments expect downloaded data!

### 0. Malignant signature extraction: CRC and ESCC
This step can be skipped if the anntations have been downloaded successfully from the Drive, see [Data Setup](#data-setup).

Run:
```bash
cd data/sh_files
bash run_dgex_cancer_sigs_with_pseudobulks.sh
bash run_dgex_cancer_sigs_on_individual_samples.sh
bash run_dgex_non_rel_genes_with_pseudobulks.sh
```

The scripts are supposed to load the preprocessed data, compute the malignant cell-specific signatures and store them in the `[Project location]/data/annotations` folder. 

### 1. Comparison of score equality between Python and R implementations
To compare the euqality of implementations in Python and R (**Figure 1b**), we follow the following steps:

1. First, load the data and convert it to single-cell experiments in R, i.e., running first part of `notebooks/construction_scoring_methods/compare_python_and_R_versions_of_scoring_methods.ipynb`
2. Adapt the missing paths in the `R` script `notebooks/construction_scoring_methods/scoring_crc_escc_luad_w_ucell_jasmine_seurat_ans.R` and run it
3. Run the rest of the notebook in 1. It will create the subplot (**Figure 1b**)

### 2. Optimal control gene selection
To generate **Figure 1c**, first run the following experiment script

```bash
cd experiments
bash run_all_experiments_on_dataset.sh
```

For the visualization:
1. **Figure 1c. and S3**: Run `notebooks/control_genes_selection_experiments/figures_contol_bias_experiment.ipynb` 
   > [!NOTE] The `test_folder` variable has to be set to a results folder created by the previous script, e.g., `[path to project]/experiments/control_genes_selection/mean_var_per_gene_scores/B memory kappa`
2. **Figure S1**: Run `notebooks/control_genes_selection_experiments/Control_gene_selection_comparison_last_expression_bin.ipynb`


### 3. Evaluating scoring method robustness to batch effects: individual vs. joint sample scoring
To get **Figure 1d**. we please run, 
```bash
cd experiments/data_composition_experiments

bash run_data_comp_exp.sh crc
bash run_data_comp_exp.sh escc
```

It will store the plots at `[project location]/experiments/data_composition_experiments/[crc | escc]/mean_norm/dgex_on_pseudobulk/strip_plots`

For **Figure S5**. run `notebooks/data_composition_experiments/variance_decrease.ipynb` 

### 4. Evaluating scoring method robustness to scoring small signatures
To recreate **Figure S4b and S6a**, first run the scoring experiments: 

```bash
cd  experiments/signature_lengths_experiments
bash run_sig_length_exp.sh escc 100 
bash run_sig_length_exp.sh crc 150 
```

Then use notebook `notebooks/signature_lengths_experiments/result_heatmaps_figures.ipynb` to create the figures.


### 5. Evaluating scoring method robustness to noise in gene expression signatures
To recreate **Figure S4c and S6b**, first run the scoring experiments: 

```bash
cd  experiments/signature_noise_addition_experiments
bash run_sig_noise_exp.sh escc 100 
bash run_sig_noise_exp.sh crc 100 
```

Then use notebook `notebooks/signature_noise_addition_experiments/sig_noise_experiment_figures.ipynb` to create the figures.


### 6. Evaluating score range comparability between scoring methods for cell state annotation
To create the subplots of **Figure 2** as well as the **Figures S7-S10**, run the following experiments scripts:

```bash
cd experiments/comparable_score_ranges
bash run_comp_range_exp.sh
```

To create the all figures run: `notebooks/comparable_score_ranges/create_plot.ipynb`. 

### 7. Runtime experiments
To rerun the runtime experiment (**Figure S2**) please run the notebook `/experiments/runtime/runtime_comparison.ipynb`.

### Case study: EMT signal decoupling
See `notebooks/EMT_scoring_experiments/` for the main case study demonstrating decoupling of EMT signals in stromal and cancer cells. It contains a README.md file with the steps to reproduce the experiments.


## 📧 Correspondance

For questions, or any encountered issues do not hesitate to contact: [Laure Ciernik](mailto:laure.ciernik@gmail.com) and [Agnieszka Kraft](mailto:agnieszka.kraft@inf.ethz.ch)


## 📜 Citation

If you use this code or the ANS method in your research, please cite:

```bibtex
@article {ans_gene_signature_scoring,
	author = {Ciernik, Laure and Kraft, Agnieszka and Barkmann, Florian and Yates, Josephine and Boeva, Valentina},
	title = {Robust and efficient annotation of cell states through gene signature scoring},
	year = {2025},
	doi = {10.1101/2023.09.20.558114},
	journal = {bioRxiv}
}
```

