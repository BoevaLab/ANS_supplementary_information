import argparse
import json
import os
import random
import sys
import warnings
from datetime import datetime

import pandas as pd
import scanpy as sc
from signaturescoring import score_signature
from signaturescoring.utils.utils import get_mean_and_variance_gene_expression
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

sys.path.append('../..')
from data.load_data import load_datasets
from data.constants import CANCER_DATASETS, NORM_METHODS, METHOD_WO_MEAN, DGEX_GRANULARITY, BASE_PATH_EXPERIMENTS
from experiments.experiment_utils import AttributeDict, get_scoring_method_params, get_malignant_signature

sc.settings.verbosity = 2


def generate_storing_path(base_path, dataset):
    st_path_dec = os.path.join(base_path, dataset, 'decreasing_log2fc', 'AUCROCS')

    if not os.path.exists(st_path_dec):
        os.makedirs(st_path_dec)
        sc.logging.info(f'> Created new directory with path {st_path_dec}')

    return st_path_dec


def score_genes_and_evaluate(adata, gene_list, df_mean_var, sc_method_long, sc_method, scm_params, col_sid='sample_id'):
    if sc_method in METHOD_WO_MEAN:
        score_signature(
            method=sc_method,
            adata=adata,
            gene_list=gene_list,
            **scm_params
        )
    else:
        score_signature(
            method=sc_method,
            adata=adata,
            gene_list=gene_list,
            df_mean_var=df_mean_var,
            **scm_params
        )
    curr_scores = adata.obs[scm_params['score_name']].copy()
    aucs = []
    precision, recall, thresholds = precision_recall_curve(adata.obs.malignant_key, curr_scores, pos_label='malignant')
    # calculate precision-recall AUC
    res_auc = auc(recall, precision)
    
    aucs.append((len(gene_list),
                 1 - roc_auc_score(adata.obs.malignant_key, curr_scores), 
                 res_auc))
        
    return pd.DataFrame(aucs, columns=['signature_length',f'AUCROC_{sc_method_long}', f'AUCPR_{sc_method_long}'])

def run_experiment_decreasing_dgex(adata, gene_list, df_mean_var, sc_method_long, sc_method, scm_params):
    results = []
    max_sig_length = len(gene_list)
    sc.settings.verbosity = 0
    for curr_sig_len in range(1, max_sig_length + 1):
        curr_gene_list = gene_list[0:curr_sig_len]
        aucs = score_genes_and_evaluate(adata, curr_gene_list, df_mean_var, sc_method_long, sc_method, scm_params)
        results.append(aucs)
    sc.settings.verbosity = 2
    results = pd.concat(results, axis=0)
#     results = pd.pivot(results, columns='sample_id', index='signature_length')
    return results


def only_allowed_genes(df_mean_var, gene_list):
    ctrl_size = 100
    gene_means = df_mean_var['mean'].copy()
    sorted_gene_means = gene_means.sort_values()
    sc.logging.info(type(sorted_gene_means.iloc[-(ctrl_size//2):]))
    not_allowed = set(sorted_gene_means.iloc[-(ctrl_size//2):].index.tolist() + sorted_gene_means.iloc[:(ctrl_size//2)].index.tolist())
    gene_list = [x for x in gene_list if not (x in not_allowed)]
    return gene_list


def main(config):
    sc.logging.info(f'Generate storing path.')
#     st_path_dec, st_path_rand = generate_storing_path(config.base_storing_path, config.dataset)
    st_path_dec = generate_storing_path(config.base_storing_path, config.dataset)

    config['st_path_dec'] = st_path_dec
#     config['st_path_rand'] = st_path_rand
    sc.logging.info(f'-----------------------------------------------------------------')
#     sc.logging.info(f'Got the storing paths {(st_path_dec, st_path_rand)}.\nStoring experiment configuration ...')
    sc.logging.info(f'Got the storing paths {st_path_dec}.\nStoring experiment configuration ...')
    start = datetime.now()
    with open(os.path.join(st_path_dec, f'config_{config.scoring_method}.txt'), "w") as fp:
        json.dump(config, fp, indent=4)  # encode dict into JSON
    # with open(os.path.join(st_path_rand, f'config_{config.scoring_method}.txt'), "w") as fp:
    #     json.dump(config, fp, indent=4)  # encode dict into JSON
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Get {config.scoring_method} parameters.')
    start = datetime.now()
    (sc_method, scm_params) = get_scoring_method_params(config.scoring_method)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Load {config.dataset} dataset.')
    start = datetime.now()
    adata = load_datasets(config.dataset, norm_method=config.norm_method)
    if 'log1p' in adata.uns_keys():
        adata.uns['log1p']['base'] = None
    else:
        adata.uns['log1p'] = {'base': None}
        
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info('Get mean and variance gene expression')
    start = datetime.now()
    df_mean_var = get_mean_and_variance_gene_expression(adata, estim_var=False)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Running experiment with signature growing in decreasing logFC order.')
    start = datetime.now()

    gene_list = get_malignant_signature(config.dataset,
                                        config.norm_method,
                                        False,
                                        config.dge_on_all,
                                        config.intersect_pct,
                                        config.min_log2fc,
                                        config.adj_pval,
                                        config.ranked_means,
                                        config.sort_log2FC_by,
                                        config.max_sig_length,
                                        most_dge=True)
    sc.logging.info(f'> len(gene_list)={len(gene_list)} before removing not allowed genes')
    gene_list = only_allowed_genes(df_mean_var, gene_list)
    sc.logging.info(f'> len(gene_list)={len(gene_list)} after removing not allowed genes')
    results = run_experiment_decreasing_dgex(adata, gene_list, df_mean_var, config.scoring_method, sc_method, scm_params)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Storing results with signature in decreasing logFC order')
    start = datetime.now()
    fn = os.path.join(st_path_dec, f'AUC_decreasing_log2fc_{config.scoring_method}.csv')
    sc.logging.info(f'Storing simulation AUCROC results ({fn}).')
    results.to_csv(fn)
    sc.logging.info(f'> Duration {datetime.now() - start}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run signature length experiment.')

    # DATASET CONFIGURATION
    parser.add_argument('--dataset', choices=CANCER_DATASETS, default='crc', help='Dataset to do the experiment on.')
    parser.add_argument('--norm_method', choices=NORM_METHODS, default='mean', help='Which normalization method used.')
    parser.add_argument('--dge_on_all', default='pseudobulk', choices=DGEX_GRANULARITY,
                        help='Use signature for malignant cells computed on all samples. Default use signature '
                             'computed on each sample individually and take genes that appear in intersect_pct % of '
                             'the samples.')
    parser.add_argument('--ranked_means', action='store_true', help='Use dgex genes sorted by mean rank log2FC. Only '
                                                                    'used if dge_on_all=individual.')
    parser.add_argument('--intersect_pct', default=0.9, type=float, help="Defines the percentage of samples a gene "
                                                                         "needs to appear in the per sample DGEX to "
                                                                         "be considered as overall DGEX. Only "
                                                                         "important for dge_on_all=individual with "
                                                                         "ranked_means=False.")
    parser.add_argument('--min_log2fc', default=2, type=float, help='Log2FC threshold for DGEX gene selection.')
    parser.add_argument('--adj_pval', default=0.01, type=float, help='Adjusted p-value threshold for DGEX gene '
                                                                     'selection.')
    parser.add_argument('--sort_log2FC_by', default='median_log2FC', choices=['mean_log2FC', 'median_log2FC'],
                        help="Defines the order in which to sort the DGEX genes found in several samples, "
                             "when dge_on_all=False.")

    # EXPERIMENT CONFIGURATION
    parser.add_argument('--scoring_method', choices=['adjusted_neighborhood_scoring',
                                                     'seurat_scoring',
                                                     'seurat_ag_scoring',
                                                     'seurat_lvg_scoring',
                                                     'scanpy_scoring',
                                                     'jasmine_scoring_lh',
                                                     'jasmine_scoring_or',
                                                     'ucell_scoring'],
                        default='adjusted_neighborhood_scoring', help='Do experiment with a specific scoring method.')

    parser.add_argument('--max_sig_length', default=250, type=int,
                        help='Maximum length of signature for malignant cells.')

    # STORING CONFIGURATION
    parser.add_argument('--base_storing_path',
                        default=os.path.join(BASE_PATH_EXPERIMENTS,'signature_lengths_experiments'),
                        help='Path where to store the results.')

    args = AttributeDict(vars(parser.parse_args()))

    start = datetime.now()
    sc.logging.info(f'Signature length experiment with the following configuration:\n{json.dumps(args, indent=4)}')

    main(args)

    sc.logging.info(f'Finished experiment in total {datetime.now() - start} time.')
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'-----------------------------------------------------------------')
