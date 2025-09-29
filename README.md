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

### 0. Malignant signature extraction: CRC and ESCC

### 1. Comparison of score equality between Python and R implementations

### 2. Optimal control gene selection
To generate Figure 1.c 

```bash
cd experiments
bash run_all_experiments_on_dataset.sh
```
### 3.Evaluating scoring method robustness to batch effects: individual vs. joint sample scoring

### 4. Evaluating scoring method robustness to scoring small signatures

### 5. Evaluating scoring method robustness to noise in gene expression signatures

### 6. Evaluating score range comparability between scoring methods for cell state annotaion

### Case study: EMT signal decoupling
See `notebooks/EMT_scoring_experiments/` for the main case study demonstrating decoupling of EMT signals in stromal and cancer cells.


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

