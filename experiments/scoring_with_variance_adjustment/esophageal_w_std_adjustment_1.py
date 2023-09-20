import datetime
import os
import pathlib
import pickle
import sys

import pandas as pd
import scanpy as sc

sys.path.append("../..")

from src.scoring_methods.gene_signature_scoring import score_signature
from experiments.evaluate_scoring_method import run_n_times_scoring
from src.utils.metrics import get_AUC_and_F1_performance
from src.utils.utils import get_mean_and_variance_gene_expression
from src.data.preprocess_data import preprocess
from data.old.load_esophageal_data import load_esophageal_data

############################################################
# Define global variables
############################################################
nr_sig_lengths = [20, 50, 100]
nr_sig_to_test = 20
# nr_sig_lengths = [20]
# nr_sig_to_test = 2
logfoldchange_cutoff = 3
# TODO define path to store data
base_storing_path = '../esophageal_w_std_adjustment/part1/'
pathlib.Path(base_storing_path).mkdir(parents=True, exist_ok=True)

scoring_methods = [
    {
        "scoring_method": "original_scanpy_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "score_name": "original_scanpy_scoring",
        },
    },
    {
        "scoring_method": "corrected_scanpy_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "score_name": "corrected_scanpy_scoring",
        },
    },
    {
        "scoring_method": "tirosh_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "score_name": "tirosh_scoring",
        },
    },
    {
        "scoring_method": "tirosh_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "adjust_for_gene_std": True,
            "score_name": "tirosh_scoring_std_adjust",
        },
    },
    {
        "scoring_method": "tirosh_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "adjust_for_gene_std": True,
            "adjust_for_gene_std_var_1p": True,
            "score_name": "tirosh_scoring_std_1p_adjust",
        },
    },
    {
        "scoring_method": "neighborhood_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "score_name": "neighborhood_scoring",
        },
    },
    {
        "scoring_method": "neighborhood_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "adjust_for_gene_std": True,
            "score_name": "neighborhood_scoring_std_adjust",
        },
    },
    {
        "scoring_method": "neighborhood_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "adjust_for_gene_std": True,
            "adjust_for_gene_std_var_1p": True,
            "score_name": "neighborhood_scoring_std_1p_adjust",
        },
    },
]
methods_to_pass_df_mean_var = ['adjusted_neighborhood_scoring',
                               'agcg_tirosh_scoring', 'lvcg_tirosh_scoring',
                               'neighborhood_scoring', 'tirosh_scoring']
methods_to_run_scoring_multiple_times = ['original_scanpy_scoring',
                                         'corrected_scanpy_scoring',
                                         'tirosh_scoring']
nr_simul = 50

############################################################
# Load and preprocess data
############################################################
adatas = load_esophageal_data(small=False, n_samples=54)

for key, adata in adatas.items():
    preprocess(adata,
               min_genes=700,
               min_cells=10,
               target_sum=1e4)

big_adata = sc.concat(list(adatas.values()), join='inner')

############################################################
# Compute mean, variance and estimated variance for data
############################################################
df_mean_var = get_mean_and_variance_gene_expression(
    big_adata,
    estim_var=True,
    show_plots=False,
    store_path=base_storing_path,
    store_data_prefix=''
)
############################################################
# Do differential gene expression and store matrix
############################################################
sc.tl.rank_genes_groups(
    big_adata, "healthy", method='wilcoxon', key_added='wilcoxon', tie_correct=True
)
wc = sc.get.rank_genes_groups_df(
    big_adata,
    group="unhealthy",
    key='wilcoxon',
    pval_cutoff=0.01,
    log2fc_min=0.01,
)
wc.to_csv(os.path.join(base_storing_path, 'DE_results.csv'))

############################################################
# For each signature length compute nr_sig_to_test a
############################################################
test_stats = []
for nr_sig_genes in nr_sig_lengths:
    print(f'Start scoring signatures with length {nr_sig_genes} ...')
    start = datetime.datetime.now()
    signatures = []
    for i in range(nr_sig_to_test):
        signatures.append(wc[wc['logfoldchanges'] > logfoldchange_cutoff].sample(nr_sig_genes)['names'].tolist())
    with open(os.path.join(base_storing_path, f'signtures_w_length_{str(nr_sig_genes)}'), "wb") as fp:  # Pickling
        pickle.dump(signatures, fp)
    for j, gene_list in enumerate(signatures):
        print(f'   > Start scoring {j}th signature...')
        scoring_names = []
        for sc_dict in scoring_methods:
            scoring_method = sc_dict["scoring_method"]
            sc_params = sc_dict["sc_params"].copy()
            sc_params["score_name"] = sc_params["score_name"] + f'_{str(j)}th_sig_w_len_{str(nr_sig_genes)}'
            scoring_names.append(sc_params["score_name"])
            print(f'      > Start scoring {sc_params["score_name"]} ...')
            if sc_dict["sc_params"]["score_name"] in methods_to_run_scoring_multiple_times:
                runs_df = run_n_times_scoring(
                    adata=big_adata,
                    gene_list=gene_list,
                    n_runs=nr_simul,
                    scoring_method=scoring_method,
                    scoring_params=sc_params,
                )
                big_adata.obs[sc_params["score_name"]] = runs_df.mean(axis=1)
            else:
                if scoring_method in methods_to_pass_df_mean_var:
                    score_signature(
                        method=scoring_method,
                        adata=big_adata,
                        gene_list=gene_list,
                        df_mean_var=df_mean_var,
                        **sc_params
                    )
                else:
                    score_signature(
                        method=scoring_method,
                        adata=big_adata,
                        gene_list=gene_list,
                        **sc_params
                    )
            print(f'      > ... Finished scoring {sc_params["score_name"]}')
        print(f'   > ... Finished scoring {j}th signature')
        print(f'   > Start measuring performance for {j}th signature scoring ...')
        test_stat = pd.DataFrame()
        for group in big_adata.obs.groupby(by='sample_id'):
            test_stat = get_AUC_and_F1_performance(
                big_adata[group[1].index, :],
                scoring_names,
                label_col="healthy",
                label_of_interest="unhealthy",
                old_df=test_stat,
                sample_id=group[0],
            )
        test_stat['sig_length'] = nr_sig_genes
        test_stat['nr_signature'] = j
        test_stats.append(test_stat)
        print(f'   > ... Finished measuring performance for {j}th signature scoring')
    end = datetime.datetime.now()
    delta = end - start
    print(f'... Finished scoring signatures with length {nr_sig_genes} (elapsed time {str(delta)})')

print('Storing big_adata.obs and test_statistics!')
big_adata.obs.to_csv(os.path.join(base_storing_path, 'scores.csv'))
all_test_stats = pd.concat(test_stats)
all_test_stats.to_csv(os.path.join(base_storing_path, 'all_test_stats.csv'))
