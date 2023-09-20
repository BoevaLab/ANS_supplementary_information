import argparse
import json
import os
import random
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
from datetime import datetime

import pandas as pd
import scanpy as sc

from signaturescoring import score_signature
from signaturescoring.utils.utils import get_mean_and_variance_gene_expression
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

sys.path.append('../..')
from data.load_data import load_datasets, load_non_rel_genes_for_mal_cells
from data.constants import BASE_PATH_EXPERIMENTS, CANCER_DATASETS, NORM_METHODS, DGEX_GRANULARITY, METHOD_WO_MEAN
from experiments.experiment_utils import AttributeDict, get_scoring_method_params, get_malignant_signature

sc.settings.verbosity = 2


def generate_storing_path(base_path, dataset, max_sig_len=100):
    st_path = os.path.join(base_path, dataset, 'AUCROCS', f'sig_len_{max_sig_len}')

    if not os.path.exists(st_path):
        os.makedirs(st_path)
        sc.logging.info(f'Created new directory with path {st_path}')

    return st_path


def get_remaining_genes(dataset, adata, gene_list, norm_method, sample_based, dge_on_all, max_abs_log2fc=None,
                        min_pval=None):
    if max_abs_log2fc is None and min_pval is None:
        all_genes = set(adata.var_names)
        return list(all_genes.difference(gene_list))

    wc = load_non_rel_genes_for_mal_cells(dataset,
                                          norm_method,
                                          sample_based,
                                          dge_on_all)
    if max_abs_log2fc is None:
        filter_g = wc.padj > min_pval
    elif min_pval is None:
        filter_g = wc.log2FoldChange.abs() <= max_abs_log2fc
    else:
        filter_g = (wc.log2FoldChange.abs() <= max_abs_log2fc) & (wc.padj > min_pval)

    return wc[filter_g].genes.tolist()


def only_allowed_genes(df_mean_var, gene_list):
    ctrl_size = 100
    gene_means = df_mean_var['mean'].copy()
    sorted_gene_means = gene_means.sort_values()
    # sc.logging.info(type(sorted_gene_means.iloc[-(ctrl_size//2):]))
    not_allowed = set(sorted_gene_means.iloc[-(ctrl_size // 2):].index.tolist() + sorted_gene_means.iloc[
                                                                                  :(ctrl_size // 2)].index.tolist())
    gene_list = [x for x in gene_list if not (x in not_allowed)]
    return gene_list


def run_one_experiment(run_id, sc_method, scm_params, adata, df_mean_var, gene_list, remaining_genes, step=0.05):
    aucrocs = []
    for i in range(0, len(gene_list) + 1, round(len(gene_list) * step)):
        nr_noise_genes = i
        nr_orig_sig_genes = len(gene_list) - i
        curr_sig = random.sample(gene_list, nr_orig_sig_genes) + random.sample(remaining_genes, nr_noise_genes)

        if sc_method in METHOD_WO_MEAN:
            score_signature(adata=adata, gene_list=curr_sig, method=sc_method, **scm_params)
        else:
            score_signature(adata=adata, gene_list=curr_sig, method=sc_method, df_mean_var=df_mean_var, **scm_params)

        curr_scores = adata.obs[scm_params['score_name']].copy()
        # AUCROC
        aucroc = 1 - roc_auc_score(adata.obs.malignant_key, curr_scores)
        # AUCPR
        precision, recall, thresholds = precision_recall_curve(adata.obs.malignant_key, curr_scores,
                                                               pos_label='malignant')
        res_auc = auc(recall, precision)
        # append results
        aucrocs.append((nr_orig_sig_genes, aucroc, res_auc))

    if i < len(gene_list):
        # score only with random genes
        nr_noise_genes = len(gene_list)
        nr_orig_sig_genes = 0
        curr_sig = random.sample(gene_list, nr_orig_sig_genes) + random.sample(remaining_genes, nr_noise_genes)

        if sc_method in METHOD_WO_MEAN:
            score_signature(adata=adata, gene_list=curr_sig, method=sc_method, **scm_params)
        else:
            score_signature(adata=adata, gene_list=curr_sig, method=sc_method, df_mean_var=df_mean_var, **scm_params)
        curr_scores = adata.obs[scm_params['score_name']].copy()
        # AUCROC
        aucroc = 1 - roc_auc_score(adata.obs.malignant_key, curr_scores)
        # AUCPR
        precision, recall, thresholds = precision_recall_curve(adata.obs.malignant_key, curr_scores,
                                                               pos_label='malignant')
        res_auc = auc(recall, precision)
        # append results
        aucrocs.append((nr_orig_sig_genes, aucroc, res_auc))

    aucrocs = pd.DataFrame(aucrocs, columns=['purity', f'AUCROC_{run_id}', f'AUCPR_{run_id}'])
    aucrocs = aucrocs.set_index('purity')
    return aucrocs


def main(config):
    sc.logging.info(f'Starting experiment ...')
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Get {config.scoring_method} parameters.')
    start = datetime.now()
    (sc_method, scm_params) = get_scoring_method_params(config.scoring_method)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Generate storing path.')
    start = datetime.now()
    storing_path = generate_storing_path(config.base_storing_path,
                                         config.dataset,
                                         config.signature_length)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    config['storing_path'] = storing_path

    sc.logging.info(f'-----------------------------------------------------------------')
    if config.max_abs_log2fc is None and config.min_pval is None:
        suffix = 'noise_genes_all_remaining'
    elif config.max_abs_log2fc is None and config.min_pval is not None:
        suffix = f'noise_genes_min_pval_{config.min_pval}'
    elif config.max_abs_log2fc is not None and config.min_pval is None:
        suffix = f'noise_genes_max_abs_log2fc_{config.max_abs_log2fc}'
    else:
        suffix = f'noise_genes_max_abs_log2fc_{config.max_abs_log2fc}_min_pval_{config.min_pval}'
    sc.logging.info(f'> File suffix {suffix}.')
    config['noise'] = suffix

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Got the storing paths {storing_path}.\nStoring experiment configuration ...')
    start = datetime.now()
    with open(os.path.join(storing_path, f'config_{config.scoring_method}_{suffix}.txt'), "w") as fp:
        json.dump(config, fp, indent=4)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Load {config.dataset} dataset.')
    start = datetime.now()
    adata = load_datasets(config.dataset,
                          norm_method=config.norm_method)
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
    sc.logging.info(f'Load differentially expressed genes for the dataset.')
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
                                        config.signature_length,
                                        most_dge=True)

    remaining_genes = get_remaining_genes(config.dataset,
                                          adata,
                                          gene_list,
                                          config.norm_method,
                                          False,
                                          config.dge_on_all,
                                          config.max_abs_log2fc,
                                          config.min_pval)
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Remove genes at the borders.')
    start = datetime.now()

    sc.logging.info(f'> len(gene_list)={len(gene_list)}, len(remaining_genes)={len(remaining_genes)} '
                    f'before removing not allowed genes')
    gene_list = only_allowed_genes(df_mean_var, gene_list)
    remaining_genes = only_allowed_genes(df_mean_var, remaining_genes)
    sc.logging.info(f'> len(gene_list)={len(gene_list)}, len(remaining_genes)={len(remaining_genes)} '
                    f'after removing not allowed genes')
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Running experiment for {config.nr_sims} simulations.')
    start = datetime.now()
    sim_res = []
    for i in range(config.nr_sims):
        sc.logging.info(f'> Starting {i}/{config.nr_sims} (passed time {datetime.now() - start}).')
        sc.settings.verbosity = 0
        res_one_run = run_one_experiment(i, sc_method, scm_params, adata, df_mean_var, gene_list, remaining_genes)
        sc.settings.verbosity = 2
        sim_res.append(res_one_run)

    sim_res = pd.concat(sim_res, axis=1)

    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    fn = os.path.join(storing_path, f'AUCROC_{config.nr_sims}sims_{config.scoring_method}_{suffix}.csv')
    sc.logging.info(f'Storing simulation AUCROC results ({config.nr_sims}) at {fn}.')
    sim_res.to_csv(fn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run signature noise-robustness experiment.')

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
    parser.add_argument('--adj_pval', default=0.01, type=float,
                        help='Adjusted p-value threshold for DGEX gene selection.')
    parser.add_argument('--sort_log2FC_by', default='median_log2FC', choices=['mean_log2FC', 'median_log2FC'],
                        help="Defines the order in which to sort the DGEX genes found in several samples, when "
                             "dge_on_all=False.")

    parser.add_argument('--max_abs_log2fc', default=None, type=float, help='Log2FC threshold for non-relevant genes.')
    parser.add_argument('--min_pval', default=None, type=float,
                        help='Adjusted p-value threshold for non-relevant genes.')

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
    parser.add_argument('--signature_length',
                        default=100, type=int, help='Length of signature for malignant cells.')
    parser.add_argument('--nr_sims',
                        default=20, type=int, help='Dataset to do the experiment on.')

    # STORING CONFIGURATION
    parser.add_argument('--base_storing_path',
                        default=os.path.join(BASE_PATH_EXPERIMENTS,'signature_noise_addition_experiments'),
                        help='Path where store the results.')

    args = AttributeDict(vars(parser.parse_args()))

    start = datetime.now()
    sc.logging.info(
        f'Signature signal-to-noise experiment with the following configuration:\n{json.dumps(args, indent=4)}')

    main(args)

    sc.logging.info(f'Finished experiment in total {datetime.now() - start} time.')
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'-----------------------------------------------------------------')
