import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from tqdm import tqdm

sys.path.append('..')
sys.path.append('../..')

from load_data import load_datasets
from constants import BASE_PATH_DGEX_CANCER
from experiments.experiment_utils import AttributeDict

# Global settings
sc.settings.verbosity = 2
plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 14})


def get_storing_path(base_storing_path, dataset, norm_method, min_log2fc, pval):
    storing_path = os.path.join(base_storing_path,
                                dataset,
                                f'{norm_method}_norm',
                                f'dgex_on_each_sid',
                                f'min_log2fc_{min_log2fc}_pval_{pval}')
    if not os.path.exists(storing_path):
        os.makedirs(storing_path)
    sc.logging.info(f'> Selected storing path {storing_path}')
    return storing_path


def get_dgex_genes(adata, logfc_min=2, pval_max=0.01):
    curr_adata = adata.copy()
    sc.tl.rank_genes_groups(curr_adata, 'malignant_key', method='wilcoxon', key_added='wilcoxon', tie_correct=True)
    # get all the genes and only select genes afterward
    wc = sc.get.rank_genes_groups_df(curr_adata, group='malignant', key='wilcoxon')
    gex_genes = wc[(wc.logfoldchanges > logfc_min) & (wc.pvals_adj < pval_max)]
    gex_genes = gex_genes.sort_values(by='logfoldchanges', ascending=False).reset_index(drop=True)
    return gex_genes


def get_per_sample_dgex_genes(adata, logfc_min=2, pval_max=0.01, col_sid='sample_id'):
    sc.settings.verbosity = 0
    adatas = {}
    for group in adata.obs.groupby(col_sid):
        adatas[group[0]] = adata[group[1].index,].copy()
        sc.pp.filter_genes(adatas[group[0]], min_cells=10)

    list_dges = []
    for sid, curr_adata in tqdm(adatas.items()):
        curr_genes = get_dgex_genes(curr_adata, logfc_min, pval_max)[['names', 'logfoldchanges']].copy()
        curr_genes = curr_genes.set_index('names')
        list_dges.append(curr_genes)
    sc.settings.verbosity = 2
    return list_dges


def get_genes_dgex_genes_in_pct_samples(list_dges):
    logfc_per_sample_and_gene = pd.concat(list_dges, axis=1, join='outer')
    logfc_per_sample_and_gene.fillna(0, inplace=True)
    ranked_logfc_per_sample_and_gene = logfc_per_sample_and_gene.rank(ascending=False, axis=0)

    mean_ranked_logfc_per_sample_and_gene = ranked_logfc_per_sample_and_gene.mean(axis=1)
    mean_ranked_logfc_per_sample_and_gene.name = 'mean_ranked_log2FC'
    sorted_mean_ranked_logfc_per_sample_and_gene = mean_ranked_logfc_per_sample_and_gene.sort_values(ascending=True)
    return sorted_mean_ranked_logfc_per_sample_and_gene.reset_index()



def main(config):
    # get storing path
    sc.logging.info(f'Get storing path.')
    storing_path = get_storing_path(config.base_storing_path, config.dataset, config.norm_method, config.min_log2fc,
                                    config.max_padj)
    config['storing_path'] = storing_path
    sc.logging.info(f'Storing experiment configuration...')
    with open(os.path.join(storing_path, 'config.txt'), "w") as fp:
        json.dump(config, fp, indent=4)  # encode dict into JSON

    # load data
    sc.logging.info(f'Load dataset.')
    adata = load_datasets(config.dataset, preprocessed=True, norm_method=config.norm_method)
    adata.uns['log1p']['base'] = None

    adata.layers['normalized'] = adata.X.copy()

    sc.logging.info(f'Compute DGEX per sample.')
    list_dges = get_per_sample_dgex_genes(adata, logfc_min=config.min_log2fc, pval_max=config.max_padj,
                                          col_sid='sample_id')

    sc.logging.info(f'Get genes sorted by mean rank of log2FC.')
    dgex_genes = get_genes_dgex_genes_in_pct_samples(list_dges)

    dgex_genes.to_csv(os.path.join(storing_path, f'dgex_genes_ranked_means.csv'), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run dataset composition experiment.')

    # DATASET CONFIGURATION
    parser.add_argument('--dataset', choices=['crc', 'escc', 'luad'],
                        default='crc', help='indicate dataset on which to create signature.')
    parser.add_argument('--norm_method', choices=['mean', 'median', 'CP10k'],
                        default='mean', help='Indicate normalization method used in the preprocessing.')
    # DGEX CONFIGURATION
    parser.add_argument('--min_log2fc', default=2, type=float, help='Minimum log2FC for DGEX.')
    parser.add_argument('--max_padj', default=0.01, type=float, help='Maximum adjusted P-value for DGEX.')

    # STORING CONFIGURATION
    parser.add_argument('--base_storing_path',
                        default=BASE_PATH_DGEX_CANCER,
                        help='Path where to store the results.')

    args = AttributeDict(vars(parser.parse_args()))

    sc.logging.info(f'Creating malignant signatures with pseudobulks with the following '
                    f'configuration:\n{json.dumps(args, indent=4)}')

    main(args)
