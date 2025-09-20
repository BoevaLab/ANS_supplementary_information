# Robust and efficient annotation of cell states through gene signature scoring
This repository accompanies the work:
Laure Ciernik*, Agnieszka Kraft*, Florian Barkmann, Joséphine Yates, and Valentina Boeva, “Robust and efficient annotation of cell states through gene signature scoring”. doi: [https://doi.org/10.1101/2023.09.20.558114](https://doi.org/10.1101/2023.09.20.558114)

*Equal contribution

It contains the Jupyter notebooks and scripts used for the experiments and visualizations. The code presented is shared for 
reproducibility purposes and is **not organized as a package**. The implementation of the novel signature scoring method ANS and others used in the article can be found and installed from the following [packaged GitHub repository](https://github.com/BoevaLab/ANS_signature_scoring). For further information, contact 
[laure.ciernik@gmail.com](mailto:laure.ciernik@gmail.com). 

### Disclaimer
After downloading the data (see below), fill in the path placeholders in `data/constants.py` as indicated by the `TODO` string. Note that the code expects the following structure: the path to the project folder should be assigned to `BASE_PATH_DRIVE`, with subfolders `data` (containing the downloaded data) and `experiments`. Placeholders for file paths have also been inserted in some jupyter notebooks, `.R` and `.sh` files. These placeholders are indicated by `TODO` strings and must be replaced for proper execution. In case of problems, do not hesitate to contact
[laure.ciernik@gmail.com](mailto:laure.ciernik@gmail.com).

## Content
The repository is structured as follows:
```
├── README.md
├── data
│   ├── ...
│   ├── constants.py
│   ├── load_data.py
│   └── run_preprocessing_cancer.sh
├── experiments
│   ├── ...
│   ├── comparable_score_ranges
│   ├── control_bias
│   ├── data_composition_experiments
│   ├── runtime
│   ├── signature_lengths_experiments
│   ├── signature_noise_addition_experiments
│   ├── run_all_experiments_on_dataset.sh
│   └── run_pbmc_experiments.sh
├── notebooks
│   ├── EMT_scoring_experiments
│   ├── comparable_score_ranges
│   ├── construction_scoring_methods
│   ├── control_genes_selection_experiments
│   ├── correlation_scores_with_TRC_and_MTP_experiments
│   ├── data_composition_experiments
│   ├── signature_lengths_experiments
│   └── signature_noise_addition_experiments
├── environment.yml
├── environment_with_build_hash.yml
└── requirements.txt
```
- The folder `data` contains helper files for preprocessing and loading desired datasets and the 
  files for DGEX to get signatures for the malignant phenotype. 
  **NOTE**: Make sure to adapt the paths in `constants.py` to the location where the data is stored. 
  See "Data Availability" below.
- The folders `experiments` and `notebooks` contain the code for executed experiments.
  - Benchmarking experiments can be run with the `experiments/run_all_cancer_experiments_on_dataset.sh` and `experiments/run_pbmc_experiments.sh` bash scripts and the jupyter notebook in the `runtime` folder. 
  - Details of the benchmarking experiments can be found in the folders of the corresponding experiment.
  - `EMT_scoring_experiments` code for the case study: decoupling the EMT signal in stromal and cancer cells. 
  
## Environment setup 
We recommend using a new Python virtual environment. We used Python 3.9.18 and a 
[miniconda](https://docs.conda.io/en/latest/miniconda.html) environment that can be created as follows:
```
conda create -n scoring_env python=3.9.18
conda activate scoring_env
pip install -r requirements.txt
```
Alternatively, you can use the `environment.yml` file to create the environment with conda (with and without build hashes).

**Note I**: There are also `environment.yml` files that can be used to create the environment with conda (with and without build hashes).

**Note II**: For running the experiments with `R`, see the `session_info.txt` file in `notebooks/construction_scoring_methods` for the used R packages.


## Data Availability
We make the preprocessed datasets (CRC, ESCC, LUAD (Xing), LUAD (Kim), breast cancer, skin cancer, ovarian cancer, and PBMC) and the used signatures available for download. There are two options to download the data:
1. Use the `data/download_preprocessed_datasets.py` script to download the data from the Google Drive folder.
2. Download the data manually from the Google Drive folder, 
[here](https://drive.google.com/drive/folders/10L2gqapJbyOn_MbrZRHQG--n0Xj7wIyg?usp=sharing). 

The Google Drive folder also contains the CanSig preprocessed datasets. In case you want to use them and get to the final preprocessed datasets, please run:
```
cd ./data/sh_files
bash run_preprocessing_cancer.sh
```


## Rerunning experiments
After downloading the data (see below) and correctly modifying the paths in `data/constants.py`, we can run the experiments. The figures for the experiments are created in the `notebooks` folder. **Note:** Manual correction of the storing paths in the jupyter notebooks might be required.  
- Preprocessing (This step can be skipped if `preprocessed` data and signatures have been downloaded): 
    - Make sure to download all the datasets of the `cansig_preprocessed` folder (see "Data Availability" above). 
    - Run the `data/run_preprocessing_cancer.sh` script for a desired dataset (CRC, ESCC, LUAD_xing, LUAD, BRCA, HGSOC, sCC). 
    - Then, get the malignant signatures by running `data/run_dgex_cancer_sigs_with_pseudobulks.sh`. 
- **TODO adapt this section as it is outdated!!!!!** Benchmark experiments (**NOTE: Change the storing paths in the bash files**):
    - The optimal control gene selection/ induced control bias experiment can be rerun by running `experiments/run_pbmc_experiments.sh` bash script. **NOTE**: Change the storing paths in the bash files of the comparable score range experiments, i.e.,`experiments/comparable_score_ranges/run_exp_[...].sh`
    - The "data composition," "signature length," and "noise robustness experiments" can be rerun by running `experiments/run_all_cancer_experiments_on_dataset.sh` bash script with a desired preprocessed dataset (CRC, ESCC, LUAD_xing).
    - The `runtime` experiment can be reproduced by running the jupyter notebook in the `experiments/runtime` folder.
    - The experiment analysing the score ranges of the different scoring methods can be rerun by running the script `experiments/run_comp_range_exp.sh`.
    - Run the different Jupyter notebooks in the `notebooks` folder to create the remaining visualizations. 
- EMT case study:
    - For each cancer type (CRC, ESCC, LUAD_xing, breast cancer) run the `get_cancer_emt_cells.ipynb` Jupyter notebook (`notebooks/EMT_scoring_experiments/CANCER_TYPE`)
    - Run `notebooks/EMT_scoring_experiments/ESCC/find_cancer_emt_signature_ESCC.ipynb`  to find the ESCC-specific cancer EMT signature 
    - Run `notebooks/EMT_scoring_experiments/LUAD_xing/find_cancer_emt_signature_LUAD.ipynb`  to find the ESCC-specific cancer EMT signature 
    - Run `notebooks/EMT_scoring_experiments/ESCC/union_ESCC_and_LUAD_specific_EMT_signature_and_refine_on_ESCC.ipynb`  to find the ESCC- and LUAD_xing-specific cancer EMT signature. 
    - Run to evaluate the new signature on CRC and breast carcinoma `notebooks/EMT_scoring_experiments/evaluation_LUNG_ESCC_signatures.ipynb`.
    - TCGA evaluation can be found in the folder `notebooks/EMT_scoring_experiments/association_with_histotypes`




## Correspondance 
First: [Laure Ciernik](mailto:laure.ciernik@gmail.com) and [Agnieszka Kraft](mailto:agnieszka.kraft@inf.ethz.ch)

Second: [Prof. Valentina Boeva](mailto:valentina.boeva@inf.ethz.ch)
